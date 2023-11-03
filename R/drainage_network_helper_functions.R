#' euclidean norm vector
#' @keywords internal
norm_vec <- function(x) sqrt(sum(x^2))

#' Finds point in distance di from point p0 in direction of point p1
#' @keywords internal
new_point <- function(p0, p1, di) {
  v = p1 - p0
  u = v / norm_vec(v)
  return (p0 + u * di)
}

#' tag nodes during creation of junctions
#' @keywords internal
tag_nodes <- function(street_vertices){

  #### tag start and end points as fixpoints ####
  street_vertices$tag <- NA
  for(i in unique(street_vertices$L1)){
    segment <- which(street_vertices$L1 == i)
    street_vertices$tag[dplyr::first(segment)] <- "fp"
    street_vertices$tag[dplyr::last(segment)] <- "fp"
  }


  #### add name for street vertices
  street_vertices$Name <- paste0("J_", seq(1:length(street_vertices$geometry)))

  # calculate distances between points
  distances <- sf::st_distance(street_vertices)

  for (i in 1:nrow(distances)){
    # select point row with distances
    check <- data.frame(dist = as.numeric(distances[i,]),
                        point = paste0("J_",seq(1:nrow(distances))), stringsAsFactors = F)

    # set zerodist to NA
    check$dist[which(check$dist == 0)] <- NA

    if(length(which(is.na(check$dist)==T)) > 1){
      # give same name for points with zero distance
      street_vertices$Name[is.na(check$dist)==T] <- check$point[is.na(check$dist)==T][1]
    }

  }


  #### tag connections ####
  street_vertices$con <- NA  # existing connections

  # check if name doubles
  for(i in 1:nrow(street_vertices)){
    n <- street_vertices$Name[i]
    if(length(which(street_vertices$Name == n)) > 1){
      street_vertices$con[i] <- T
    }else{
      street_vertices$con[i] <- F
    }
  }

  # #### tag points that need to be snapped ####
  #
  # street_vertices$mis_con <- NA  # missing connection
  #
  # # polylines need at least one fp which is connected...
  # for(i in unique(street_vertices$L1)){
  #   segment <- which(street_vertices$L1 == i)
  #   fp <- segment[which(street_vertices$tag[segment] == "fp")]
  #   if(all(street_vertices$con[fp] == F)){
  #     street_vertices$mis_con[fp] <- TRUE
  #   }
  # }

  # return the data set:
  return(street_vertices)
}

#' order tags for start and end nodes
#' @keywords internal
define_tag_order <- function(street_vertices){
  street_vertices$tag_order <- NA
  for(i in unique(street_vertices$L2)){
    segment <- which(street_vertices$L2 == i)
    street_vertices$tag_order[dplyr::first(segment)] <- "start"
    street_vertices$tag_order[dplyr::last(segment)] <- "end"
  }
  return(street_vertices)
}

#' add a new point to an existing street
#' @keywords internal
add_point_to_street <- function(point_to_add, street_vertices){
  points_in_street <- street_vertices[street_vertices$L1 == point_to_add$nearest_street,]

  # find position to add point:
  dist <- as.numeric(sf::st_distance(point_to_add, points_in_street))
  nearest_id <- which(dist == min(dist))
  second_nearest_id <- which(dist == min(dist[dist != min(dist)]))

  # order Names between which the new point should be added:
  if(nearest_id < second_nearest_id){
    first_p <- points_in_street$Name[nearest_id]
    second_p <- points_in_street$Name[second_nearest_id]
  }else{
    first_p <- points_in_street$Name[second_nearest_id]
    second_p <- points_in_street$Name[nearest_id]
  }

  # find position to add the additional point:
  id_p1 <- which(street_vertices$Name == first_p & street_vertices$L1 == point_to_add$nearest_street)
  id_p2 <- which(street_vertices$Name == second_p & street_vertices$L1 == point_to_add$nearest_street)

  # add point between nearest and second nearest vertex:
  street_vertices <- rbind(street_vertices[1:id_p1,], point_to_add, street_vertices[id_p2:nrow(street_vertices),])

  return(street_vertices)
}

#' create polylines from point data
#' @keywords internal
create_linestrings_from_vertices <- function(vertices){
  conduits <- list()
  
  layer_id <- unique(vertices$L2)

  # create line strings
  for (i in layer_id){
    # subset data based on street
    street <- vertices[vertices$L2 == i,]
    conduits[[i]] <- list()

    for(j in 1:(length(street$L2)-1)){
      # create linestrings
      conduit <- street[c(j,j+1),"L2"]
      conduit <- dplyr::summarize(dplyr::group_by(conduit, L2), .groups = "drop_last")
      # store as linestring...
      conduits[[i]][[j]] <- sf::st_cast(conduit,"LINESTRING")
    }
  }

  # unlist conduits
  conduits_L <- list()
  for( i in layer_id){
    conduits_L[[i]] <- do.call(rbind, conduits[[i]])
  }
  conduits_sf_simple <- do.call(rbind, conduits_L)

  return(conduits_sf_simple)
}

#' create polylines from point data
#' @keywords internal
create_linestrings_from_single_layer <- function(new_segment){
  new_linestring <- list()

  for(j in 1:(length(new_segment$L2)-1)){
    # create linestrings
    conduit <- new_segment[c(j,j+1),"L2"]
    conduit <- dplyr::summarize(dplyr::group_by(conduit, L2), .groups = "drop_last")
    # store as linestring...
    new_linestring[[j]] <- sf::st_cast(conduit, "LINESTRING")
  }
  new_linestring <- do.call(rbind, new_linestring)

  return(new_linestring)
}

#' Helper function for function DouglasPeuckerEpsilon copied from kmlShape Package
#' @keywords internal
shortestDistanceToLines <- function(Mx,My,Ax,Ay,Bx,By){
  aire <- abs((By-Ay)*(Mx-Ax)-(Bx-Ax)*(My-Ay))
  return(  aire / sqrt((Bx-Ax)^2 + (By-Ay)^2))
}

#' Helper function for function DouglasPeuckerEpsilon copied from kmlShape package
#' @keywords internal
findFarestPoint <- function(trajx,trajy){
  dmax <- 0
  index <- 1
  end <- length(trajx)
  
  if(end==2){
    index <- 1
    dmax <- 0
  }else{
    for(i in 2:(end-1)){
      d <- shortestDistanceToLines(Mx=trajx[i],My=trajy[i], Ax=trajx[1],Ay=trajy[1], Bx=trajx[end],By=trajy[end])
      if ( d > dmax ) {
        index <- i
        dmax <- d
      }else{}
    }
  }
  return(c(index=index,dmax=dmax))
}


#' DouglasPeuckerEpsilon function is copied from the kmlShape Package which is no longer supported
#' @keywords internal
DouglasPeuckerEpsilon <- function(trajx,trajy,epsilon,spar=NA){
  missings <- is.na(trajx)|is.na(trajy)
  if(any(missings)){
    trajx <- trajx[!missings]
    trajy <- trajy[!missings]
  }else{}
  
  if(!is.na(spar)){trajy <- smooth.spline(trajx,trajy,spar=spar)[["y"]]}else{}
  
  farestPoint <- findFarestPoint(trajx,trajy)
  index <- farestPoint["index"]
  end <- length(trajx)
  if ( farestPoint["dmax"] > epsilon ) {
    recResults1 = DouglasPeuckerEpsilon(trajx[1:index],trajy[1:index], epsilon)
    recResults2 = DouglasPeuckerEpsilon(trajx[index:end],trajy[index:end], epsilon)
    
    resultTrajx = c(recResults1$x,recResults2$x[-1])
    resultTrajy = c(recResults1$y,recResults2$y[-1])
    #        d = c(farestPoint["dmax"],recResults1$d,recResults2$d)
  } else {
    resultTrajx = c(trajx[1],trajx[end])
    resultTrajy = c(trajy[1],trajy[end])
    #        d=numeric()
  }
  return(data.frame(x=resultTrajx,y=resultTrajy))
}


#' create junctions based on street polylines
#' @keywords internal
create_junctions <- function(streets, buffer, snap_dist, epsilon, lim, junc_depth, dtm, crs_default, pre_def_junctions, pre_def_conduits){

  # plot initial data:
  graphics::plot(sf::st_geometry(streets), col = "grey", main = "street polylines")

  #### street polyline quality checks: ####
  message(" ... start quality check of street polylines")

  # extract street coordinates
  coords_streets <- as.data.frame(sf::st_coordinates(streets))

  # check layer information of streets:
  if(!(ncol(coords_streets) == 3)){
    if(length(unique(coords_streets$L1)) == nrow(streets)){
      coords_streets <- coords_streets[,c("X","Y","L1")]
    }else{
      if(length(unique(coords_streets$L2)) == nrow(streets)){
        coords_streets <- coords_streets[,c("X","Y","L2")]
        colnames(coords_streets)[3] <- "L1"
      }else{
        stop("too many polyline layers after extraction of vertices")
      }
    }
  }

  # sf of street coordinates
  street_vertices <- sf::st_as_sf(coords_streets, coords = c("X", "Y"), crs = crs_default)

  # summarize points, which belong to different layers, within a buffer :
  coords_streets$buf <- NA

  for(i in 1:nrow(street_vertices)){

    if(any(as.numeric(sf::st_distance(street_vertices[i,], street_vertices)) < buffer &
           as.numeric(sf::st_distance(street_vertices[i,], street_vertices)) != 0)){

      same <- which(as.numeric(sf::st_distance(street_vertices[i,], street_vertices)) < buffer
                    & as.numeric(sf::st_distance(street_vertices[i,], street_vertices)) != 0)
      
      # keep points that belong to different layers:
      same <- same[street_vertices$L1[same] != street_vertices$L1[i]]
      
      if(all(i < same)){
        # overwrite coordinates...
        coords_streets$X[same] <- coords_streets$X[i]
        coords_streets$Y[same] <- coords_streets$Y[i]
        coords_streets$buf[same] <- TRUE # add a tag if coordinates were adjusted
      }else{
        #... checks... can be deleted?
        dist <- sf::st_distance(street_vertices[i,], street_vertices[same,])
      }
    }

    street_vertices <- sf::st_as_sf(coords_streets, coords = c("X", "Y"), crs = crs_default)

  }

  street_vertices <- sf::st_as_sf(coords_streets, coords = c("X", "Y"), crs = crs_default)

  # plot buffered vertices:
  graphics::plot(sf::st_geometry(streets), col = "grey", main = "buffered vertices")
  graphics::plot(sf::st_geometry(street_vertices[street_vertices$buf == T,]), col = "red", add = T)

  # tag nodes for further aggregation:
  street_vertices <- tag_nodes(street_vertices)

  # plot initial fix points:
  graphics::plot(sf::st_geometry(streets), col = "grey", main = "initial fix points")
  graphics::plot(sf::st_geometry(street_vertices[street_vertices$tag == "fp",]), col = "red", add = T)

  # delete points that occur twice in one layer:
  for(i in unique(street_vertices$L1)){
    subset <- street_vertices[street_vertices$L1 == i,]
    if(any(duplicated(subset$Name))){
      # delete the point that is duplicated and has no "fp" tag:
      doubled_names <- subset$Name[duplicated(subset$Name)]

      # test if line consist only of two duplicated points:
      if(nrow(subset) <= 2){
        # delete the whole line if coordinates of both points double:
        del_node <- which(street_vertices$L1 == i)

      }else{
        for(n in doubled_names){
          if(all(is.na(subset$tag[subset$Name == n])) | all(is.na(subset$tag[subset$Name == n]) == F)){
            # if both points are NA or "fp" delete one:
            del_node <- which(street_vertices$L1 == i &
                                street_vertices$Name == n)[1]
            
          }else{
            # if only one of the doubled points is a "fp" delete the non "fp":
            del_node <- which(street_vertices$L1 == i &
                                street_vertices$Name == n &
                                is.na(street_vertices$tag) == T)
          }
        }

      }

      # delete point from data set:
      street_vertices <- street_vertices[-del_node,]

    }else{
      next
    }
  }

  # tag nodes for further aggregation:
  street_vertices <- tag_nodes(street_vertices)

  # Add vertices to near polylines:
  street_vertices$nearest_street <- NA

  for( i in which(street_vertices$con == F)){
    dist <- as.numeric(sf::st_distance(street_vertices[i,], streets))
    min_value <- min(dist[dist > 0])
    if(min_value < snap_dist){
      street_vertices$nearest_street[i] <- which(dist == min_value)[1]
    }
  }

  # plot connections between near polylines:
  graphics::plot(sf::st_geometry(streets), col = "grey", main = "missing connections")
  graphics::plot(sf::st_geometry(street_vertices[is.na(street_vertices$nearest_street) == F,]), col = "red", add = T)

  # find position to add point and add it to nearest line:
  con_name <- street_vertices$Name[which(is.na(street_vertices$nearest_street) == F)]

  for(i in con_name){
    point_to_add <- street_vertices[street_vertices$Name == i,]
    # change layer id to nearest street id:
    point_to_add$L1 <- point_to_add$nearest_street
    # delete fp tag:
    point_to_add$tag <- NA

    # add point to data:
    street_vertices <- add_point_to_street(point_to_add = point_to_add, street_vertices)

  }

  # tag nodes for further tests:
  street_vertices <- tag_nodes(street_vertices)

  # plot crossings:
  graphics::plot(sf::st_geometry(streets), col = "grey", main = "connections/crossings")
  graphics::plot(sf::st_geometry(street_vertices[street_vertices$con == T,]), col = "red", add = T)

  # add point at intersection of polylines:
  # intersect street polylines and identify crossings:
  crossings <- suppressWarnings(sf::st_intersection(streets, streets))
  crossings <- sf::st_collection_extract(crossings, "POINT")

  # plot points at intersection of polylines:
  graphics::plot(sf::st_geometry(streets), col = "grey", main = "additional nodes at line crossings")

  for(i in 1:nrow(crossings)){
    dist <- as.numeric(sf::st_distance(crossings[i,], street_vertices))
    if(all(dist > buffer)){
      # add crossings if these points do not exist in polyline coordinates...
      dist_to_streets <- as.numeric(sf::st_distance(crossings[i,], streets))
      street_L <- which(dist_to_streets < snap_dist)

      # plot additonal crossings...
      graphics::plot(sf::st_geometry(crossings[i,]), add = T, col = "red")
      graphics::plot(sf::st_geometry(streets[street_L,]), add = T, col = "red")

      # add point to crossed streets:
      for( L in street_L){
        # define crossing as point to add
        point_to_add <- crossings[i,]
        point_to_add <- as.data.frame(sf::st_coordinates(point_to_add))
        point_to_add <- unique(point_to_add)
        point_to_add <- sf::st_as_sf(point_to_add, coords = c("X", "Y"), crs = crs_default)

        point_to_add$nearest_street <- L
        point_to_add$L1 <- L
        point_to_add$tag <- NA
        point_to_add$buf <- NA
        point_to_add$con <- TRUE
        point_to_add$Name <- "crossing"

        # order point to add:
        point_to_add <- point_to_add[,c("L1", "buf", "tag", "Name", "con", "nearest_street", "geometry")]

        # add point to data:
        street_vertices <- add_point_to_street(point_to_add = point_to_add, street_vertices)
      }

    }
  }


  # double points at connections:
  while(TRUE){
    dp <- which(street_vertices$con == T & is.na(street_vertices$tag) == T)
    if(length(dp) >= 1){
      street_vertices$tag[dp[1]] <- "fp"
      street_vertices <- rbind(street_vertices[1:dp[1],], street_vertices[dp[1],], street_vertices[(dp[1]+1):nrow(street_vertices),])
    }else{
      break
    }
  }

  # update linestring IDs (L):
  street_vertices$L2 <- NA
  lfn <- 1

  while(TRUE){
    fp <- which(street_vertices$tag == "fp" & is.na(street_vertices$L2) == T)
    if(length(fp) >= 2){
      street_vertices$L2[fp[1]:fp[2]] <- lfn
      lfn <- lfn + 1
    }else{
      break
    }
  }
  
  # add crossing label to points that occur more than twice:
  for( i in unique(street_vertices$Name[street_vertices$Name != "crossing"])){
    if(length(which(street_vertices$Name == i)) > 2){
      street_vertices$Name[street_vertices$Name == i] <- "crossing"
    }
  }
  
  # update row numbers
  row.names(street_vertices) <- seq(1:nrow(street_vertices))
  
  # define start and end vertex of polyline according to order in data.frame
  street_vertices <- define_tag_order(street_vertices)
  
  # summaries polylines sharing start and end nodes (connected fixpoints, no crossings)
  vertices_to_summaries <- street_vertices$Name[street_vertices$con == "TRUE" & street_vertices$Name != "crossing" & !is.na(street_vertices$tag_order)]

  if(!is.null(vertices_to_summaries)){
    for( i in unique(vertices_to_summaries)){

      # combine only two terminating points (start/end):
      if(any(is.na(street_vertices$tag_order[street_vertices$Name ==  i]))){
        next
      }else{
        layers_to_summarize <- street_vertices$L2[street_vertices$Name == i]

        segment_to_shift <- street_vertices[street_vertices$L2 == max(layers_to_summarize),]
        segment_to_connect <- street_vertices[street_vertices$L2 == min(layers_to_summarize),]


        # check order of vertices
        if(all(street_vertices$tag_order[street_vertices$Name ==  i] == c("start", "start"))){
          # shift vertices of connecting segment (starting from second vertex) in reverse order above the connecting vertex
          combined_segments <- rbind(segment_to_shift[nrow(segment_to_shift):2,], segment_to_connect)
        }

        if(all(street_vertices$tag_order[street_vertices$Name ==  i] == c("start", "end"))){
          # shift vertices of connecting segment above the connecting vertex and keep order
          combined_segments <- rbind(segment_to_shift[1:(nrow(segment_to_shift)-1),], segment_to_connect)
        }

        if(all(street_vertices$tag_order[street_vertices$Name ==  i] == c("end", "start"))){
          # keep order but delete duplicated vertex
          combined_segments <- rbind(segment_to_connect, segment_to_shift[2:nrow(segment_to_shift),])
        }

        if(all(street_vertices$tag_order[street_vertices$Name ==  i] == c("end", "end"))){
          # shift vertices of connecting segment in reverse order below the connecting vertex
          combined_segments <- rbind(segment_to_connect, segment_to_shift[(nrow(segment_to_shift)-1):1,])
        }

        # update layer id
        combined_segments$L2 <- min(layers_to_summarize)

        # delete segments from complete dataset
        street_vertices <- street_vertices[!street_vertices$L2 %in% layers_to_summarize,]

        # rbind segments with complete dataset
        street_vertices <- rbind(street_vertices, combined_segments)

        # update tag order
        street_vertices <- define_tag_order(street_vertices)

      }
    }
  }
  
  #### Simplify shape of network: ####
  message(" ... simplify the shape of the network (Douglas-Peucker Algorithm)")

  coords_streets <- as.data.frame(sf::st_coordinates(street_vertices))
  coords_streets$L2 <- street_vertices$L2

  new_segments <- list()

  if(!is.null(pre_def_junctions)){
    new_linestrings <- list()
  }

  for(L in unique(coords_streets$L2)){
    segment <- coords_streets[coords_streets$L2 == L,]

    # Douglas-Peucker line simplification algorithm:
    simplified_segment <- DouglasPeuckerEpsilon(trajx = segment$X, trajy = segment$Y, epsilon = epsilon)

    # add point in lim distance if necessary:

    # keep first point
    result <- simplified_segment[1,]#,drop=FALSE]

    # test whether new point must be created or not
    for ( i in 1:length(simplified_segment[,1])){
      if ( i+1 < length(simplified_segment[,1]) | i+1 == length(simplified_segment[,1])){
        point <- simplified_segment[i,]
        next_point <- simplified_segment[i+1,]
        dist <- norm_vec(point[,c("x","y")] - next_point[,c("x","y")])

        if( dist/lim > 1.5){
          # move point to result if it is not added yet
          if(!all(point %in% result)){
            result <- rbind(result, point)
          }
          # create new point at lim distance as long as dist is greater than lim and move it to points
          while(TRUE){
            np <- new_point(point[,c("x","y")], next_point[,c("x","y")], lim)
            #np$L1 <- L
            result <- rbind(result, np)

            dist <- norm_vec(np[,c("x","y")] - next_point[,c("x","y")])
            if(dist/lim < 1.5){
              break
            }

            point <- np

          }
        } else {
          if(!all(point %in% result)){
            result <- rbind(result, point)
          }
        }
      }

    }

    # add last point
    result <- rbind(result, simplified_segment[nrow(simplified_segment),])

    # convert to sf:
    new_segment <- sf::st_as_sf(result, coords = c("x", "y"), crs = crs_default)

    # add layer ID...
    new_segment$L2 <- L

    # if an point object with heigths is given, transfer the heights and interpolate between them:
    if(!is.null(pre_def_junctions) & !is.null(pre_def_conduits)){

      if(!any(c("InOffset", "OutOffset") %in% colnames(pre_def_conduits))){
        pre_def_conduits$InOffset <- 0
        pre_def_conduits$OutOffset <- 0
      }

      new_segment$Bottom <- NA
      new_segment$Name <- NA
      new_segment$offset_height <- NA

      # add existing heights
      for(i in 1:nrow(new_segment)){
        id <- which(as.numeric(sf::st_distance(new_segment[i,], pre_def_junctions)) == 0)
        if(length(id) == 1){
          new_segment$Bottom[i] <- pre_def_junctions$Bottom[id]
          new_segment$Name[i] <- pre_def_junctions$Name[id]
        }
      }

      # segment length and nodes
      segment_length <- nrow(new_segment)
      first <- new_segment$Name[1]
      last <- new_segment$Name[segment_length]

      # assign names to new nodes:
      new_nodes <- length(which(is.na(new_segment$Name) == T))
      new_segment$Name[is.na(new_segment$Name)] <- paste0("J_",L,"_", 1:new_nodes)

      # pre defined conduit along which new junction heights are interpolated
      id_pre_conduit <- which((pre_def_conduits$FromNode == first & pre_def_conduits$ToNode == last) |
                                (pre_def_conduits$FromNode == last & pre_def_conduits$ToNode == first))

      # check if a predefined link exists:
      if(!length(id_pre_conduit)){
        if(any(is.na(new_segment$Bottom))){
          # if no predefined link is defined add standard Bottom height:
          new_segment$Bottom[is.na(new_segment$Bottom)] <- add_junction_heights(new_segment[is.na(new_segment$Bottom),], dtm) - junc_depth
        }

        # add links and cut link at first FromNode:
        if(new_segment$Bottom[1] > new_segment$Bottom[segment_length]){
          new_segment <- new_segment[-1,]

          # update segment length
          segment_length <- nrow(new_segment)

          if(segment_length > 1){
            # define conduits with directions
            new_linestring <- create_linestrings_from_single_layer(new_segment)
            new_linestring$FromNode <- new_segment$Name[-segment_length]
            new_linestring$ToNode <- new_segment$Name[-1]
            new_linestring$InOffset <- 0
            new_linestring$OutOffset <- 0
          }else{
            new_linestring <- NULL
            new_segment <- NULL
          }

        }else{
          new_segment <- new_segment[-segment_length,]

          # update segment length
          segment_length <- nrow(new_segment)

          if(segment_length > 1){
            # define conduits with directions
            new_linestring <- create_linestrings_from_single_layer(new_segment)
            new_linestring$FromNode <- new_segment$Name[-1]
            new_linestring$ToNode <- new_segment$Name[-segment_length]
            new_linestring$InOffset <- 0
            new_linestring$OutOffset <- 0
          }else{
            new_linestring <- NULL
            new_segment <- NULL
          }

        }

      }else{
        # if a predefined link exists interpolate heights:
        if(pre_def_conduits$FromNode[id_pre_conduit] == first){
          # assign in and outoffsets to nodes
          new_segment$offset_height[1] <- pre_def_conduits$InOffset[id_pre_conduit]
          new_segment$offset_height[segment_length] <- pre_def_conduits$OutOffset[id_pre_conduit]

          # define conduits with directions
          new_linestring <- create_linestrings_from_single_layer(new_segment)
          new_linestring$FromNode <- new_segment$Name[-segment_length]
          new_linestring$ToNode <- new_segment$Name[-1]

          # assign In- and OutOffset from pre conduits:
          new_linestring$InOffset <- 0
          new_linestring$OutOffset <- 0
          new_linestring$InOffset[1] <- pre_def_conduits$InOffset[id_pre_conduit]
          new_linestring$OutOffset[nrow(new_linestring)] <- pre_def_conduits$OutOffset[id_pre_conduit]

        }else{
          # assign in and outoffsets to nodes
          new_segment$offset_height[segment_length] <- pre_def_conduits$InOffset[id_pre_conduit]
          new_segment$offset_height[1] <- pre_def_conduits$OutOffset[id_pre_conduit]

          # define conduits with directions
          new_linestring <- create_linestrings_from_single_layer(new_segment)
          new_linestring$FromNode <- new_segment$Name[-1]
          new_linestring$ToNode <- new_segment$Name[-segment_length]

          # assign In- and OutOffset from pre conduits:
          new_linestring$InOffset <- 0
          new_linestring$OutOffset <- 0
          new_linestring$InOffset[nrow(new_linestring)] <- pre_def_conduits$InOffset[id_pre_conduit]
          new_linestring$OutOffset[1] <- pre_def_conduits$OutOffset[id_pre_conduit]

        }

        # distance weighted linear interpolation between given heights
        interpolated <- gstat::idw((Bottom + offset_height) ~ 1, new_segment[is.na(new_segment$Bottom)==F,], new_segment, idp = 1, debug.level = 0)
        new_segment$Bottom[-c(1, segment_length)] <- interpolated$var1.pred[-c(1, segment_length)]
      }

    }

    # store...
    new_segments[[L]] <- new_segment

    if(!is.null(pre_def_junctions)){
      new_linestrings[[L]] <- new_linestring
    }

  }

  street_vertices_simplified <- do.call(rbind, new_segments)

  # plot simplified network:
  graphics::plot(sf::st_geometry(streets), col = "grey", main = "nodes of simplified polylines (Douglas-Peucker Algorithm and lim distance)")
  graphics::plot(sf::st_geometry(street_vertices_simplified), col = "red", add = T)

  #### First try conduits: ####
  message(" ... define conduits along the streets without flow direction")

  if(!is.null(pre_def_junctions)){
    conduits_sf_simple <- do.call(rbind, new_linestrings)

    graphics::plot(sf::st_geometry(streets), col = "grey", main = "first try of new polylines")
    graphics::plot(sf::st_geometry(conduits_sf_simple), col = "red", add = T)
  }

  if(is.null(pre_def_junctions)){
    #conduits_sf_simple <- create_linestrings_from_vertices(vertices = street_vertices_simplified)
    
    conduits_sf_simple <-  dplyr::summarise(dplyr::group_by(street_vertices_simplified, L2), do_union = FALSE)
    conduits_sf_simple <- sf::st_cast(conduits_sf_simple, "LINESTRING")

    # plot progress...
    graphics::plot(sf::st_geometry(streets), col = "grey", main = "first try of new polylines")
    graphics::plot(sf::st_geometry(conduits_sf_simple), col = "red", add = T)

  }

  #### name junctions: ####
  message(" ... name junctions")

  junctions <- street_vertices_simplified

  if(is.null(pre_def_junctions)){

    # calculate distances between points
    distances <- sf::st_distance(junctions)

    junctions$Name <- paste0("J_", seq(1:length(junctions$geometry)))

    for (i in 1:nrow(distances)){
      # select point row with distances
      check <- data.frame(dist = as.numeric(distances[i,]),
                          point = paste0("J_",seq(1:nrow(distances))), stringsAsFactors = F)

      # set zerodist to NA
      check$dist[which(check$dist == 0)] <- NA

      if(length(which(is.na(check$dist)==T)) > 1){
        # give same name for points with zero distance
        junctions$Name[is.na(check$dist)==T] <- check$point[is.na(check$dist)==T][1]
      }

    }

  }

  if(!is.null(pre_def_junctions)){
    # keep only unique geometries:
    junctions <- dplyr::distinct(junctions, across(geometry), .keep_all = T)

  }

  #### return points of class sf ####

  if("L2" %in% colnames(junctions)){
    colnames(junctions)[colnames(junctions) == "L2"] <- "L1"
  }

  if(!is.null(pre_def_junctions)){
    if("L2" %in% colnames(conduits_sf_simple)){
      colnames(conduits_sf_simple)[colnames(conduits_sf_simple) == "L2"] <- "L1"
    }
    return(list(junctions = junctions[c("Name", "L1", "Bottom")], conduits = conduits_sf_simple))
  }else{
    return(junctions[c("Name", "L1")])
  }


}

#' add heights from dtm
#' @keywords internal
add_junction_heights <- function(junctions, dtm){
  # add heights to nodes for initial call otherwise interpolate between existing heights
  if(any(class(dtm) == "RasterLayer")){
    junctions$Top <- raster::extract(dtm, junctions)
  }else{
    if(any(class(dtm) == "sf")){
      heights <- suppressWarnings(sf::st_intersection(junctions, dtm))
      junctions$Top <- NA
      for(i in heights$Name){
        junctions$Top[junctions$Name == i] <- heights$Z[heights$Name == i]
      }
    }
  }

  if(any(is.na(junctions$Top))){
    # fill NA in top heights
    distances <- sf::st_distance(junctions)
    for(i in which(is.na(junctions$Top) == T)){
      junctions$Top[i] <- junctions$Top[which(distances[,i] == min(distances[as.numeric(distances[,i]) > 0,i]))[1]]
    }
    # if still NA fill value with mean heights:
    for(i in which(is.na(junctions$Top) == T)){
      junctions$Top[i] <- mean(junctions$Top, na.rm = T)
    }
  }

  # return junctions:
  return(junctions$Top)
}

#' create conduits along the heights of the junctions
#' @keywords internal
create_conduits_along_topography <- function(vertices, main){

  conduits <- list()
  layer_id <- unique(vertices$L1)

  # create line strings
  for (i in layer_id){
    # subset data based on street
    street <- vertices[vertices$L1 == i,]
    conduits[[i]] <- list()

    for(j in 1:(length(street$L1)-1)){
      # create linestrings
      conduit <- street[c(j,j+1),"L1"]
      conduit <- dplyr::summarise(dplyr::group_by(conduit, L1), do_union = FALSE)
      conduit <- sf::st_cast(conduit, "LINESTRING")

      # define flow directions add FromNode and ToNode
      if(street$Bottom[j] < street$Bottom[j+1]){
        conduit$ToNode <- street$Name[j]
        conduit$FromNode <- street$Name[j+1]
      }else{
        conduit$ToNode <- street$Name[j+1]
        conduit$FromNode <- street$Name[j]
      }

      # store...
      conduits[[i]][[j]] <- conduit
    }
  }

  # unlist conduits
  conduits_L <- list()
  for( i in layer_id){
    conduits_L[[i]] <- do.call(rbind, conduits[[i]])
  }
  conduits_sf <- do.call(rbind, conduits_L)

  # delete duplicated vertices
  vertices <- vertices[duplicated(vertices$geometry)==F,]

  # plot progress:
  coords <- as.data.frame(sf::st_coordinates(vertices))
  coords$Name <- vertices$Name
  # initiate plot
  graphics::plot(coords$X, coords$Y, xlab = "X", ylab = "Y", main = main)
  for(i in 1:nrow(conduits_sf)){
    p0 <- conduits_sf$FromNode[i]
    p1 <- conduits_sf$ToNode[i]
    # plot arrow:
    graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
           x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
           col = "red", length = 0.1)
  }


  # return:
  return(list(junctions = vertices, conduits_sf = conduits_sf))
}

#' add conduits to connect network to outfalls
#' @keywords internal
add_outfall_connection <- function(outfalls, junctions, conduits_sf, min_slope, junc_depth, lim){
  # add outfall conduit connection, to nearest (and higher junction)
  for( i in 1:length(outfalls$Name)){
    outfall <- outfalls[i,c("Name", "Elevation")]
    colnames(outfall) <- c("Name", "Bottom", "geometry") # change naming to group conduits later
    outfall$Top <- outfall$Bottom + junc_depth  # add column to group conduits later
    distances <- sf::st_distance(outfall, junctions)
    nearest <- sort(as.numeric(distances))[1:2]
    if(abs(nearest[1] - nearest[2]) < lim &
       junctions$Top[which(as.numeric(distances) == nearest[1])] != junctions$Top[which(as.numeric(distances) == nearest[2])]
    ){
      if(junctions$Top[which(as.numeric(distances) == nearest[1])] < junctions$Top[which(as.numeric(distances) == nearest[2])]){
        to_outfall <- junctions$Name[which(as.numeric(distances) == nearest[1])]
      }else{
        to_outfall <- junctions$Name[which(as.numeric(distances) == nearest[2])]
      }
    }else{
      to_outfall <- junctions$Name[distances == min(distances)] # nearest junction
    }

    outfall$L1 <- junctions$L1[junctions$Name == to_outfall][1]
    conduit_outfall <- rbind(outfall, junctions[junctions$Name == to_outfall,][1,])
    conduit_outfall <- dplyr::summarize(dplyr::group_by(conduit_outfall, L1), .groups = "drop_last")
    conduit_outfall <- sf::st_cast(conduit_outfall, "LINESTRING")
    conduit_outfall$FromNode <- to_outfall
    conduit_outfall$ToNode <- outfall$Name
    conduits_sf <- rbind(conduit_outfall, conduits_sf)

    # calculate elevation of outfall according to min_slope:
    outfalls$Elevation[i] <- calculate_new_bottom_height(length = as.numeric(sf::st_length(conduit_outfall)),
                                                         slope = min_slope,
                                                         FromHeight = junctions$Bottom[junctions$Name == to_outfall])

  }

  # plot progress...

  coords <- as.data.frame(sf::st_coordinates(junctions))
  coords$Name <- junctions$Name

  plot_drainage_network(coords, outfalls, conduits_sf, main = "add outfalls", col_arrow = "black")
  graphics::legend("topright", pch = 19, col = "green", c("outfalls"))

  return(list(junctions = junctions, conduits = conduits_sf, outfalls = outfalls))
}

#' network error checks
#' @keywords internal
error_checks <- function(junctions, conduits_sf, outfalls){
  junctions$tag <- NA

  for (name in unique(junctions$Name)){
    # Test 1: normal if exactly one FromNode and one ToNode
    if(length(which(conduits_sf$FromNode == name)) == 1 & length(which(conduits_sf$ToNode == name)) == 1){
      junctions$tag[junctions$Name %in% name] <- "normal"
    }
    # Test 2a: test for outfall node with exactly one ToNode and no FromNode
    if(length(which(conduits_sf$ToNode == name)) == 1 & (name %in% conduits_sf$FromNode) == F){
      junctions$tag[junctions$Name %in% name] <- "outfall_artificial"
    }
    # Test 2b: test for real outfall
    if(name %in% outfalls$Name){
      junctions$tag[junctions$Name %in% name] <- "outfall_real"
    }
    # Test 3: test for start node no ToNode one FromNode
    if((name %in% conduits_sf$ToNode) == F & length(which(conduits_sf$FromNode == name)) == 1){
      junctions$tag[junctions$Name %in% name] <- "start"
    }
    # Test 4: test for crossings or connections one FromNode more than one ToNode
    if(length(which(conduits_sf$FromNode == name)) == 1 & length(which(conduits_sf$ToNode == name)) > 1){
      junctions$tag[junctions$Name %in% name] <- "crossing"
    }

    # Test 5: sink if no FromNode and several ToNode
    if((name %in% conduits_sf$FromNode) == F & length(which(conduits_sf$ToNode == name)) > 1){
      junctions$tag[junctions$Name %in% name] <- "sink"
    }
    # Test 6: slope error if one or more ToNode and more than one FromNode
    if(length(which(conduits_sf$FromNode == name)) > 1 & length(which(conduits_sf$ToNode == name)) >= 1){
      junctions$tag[junctions$Name %in% name] <- "crossing" # slope_error
    }
    # Test 7: test for break node no ToNode several FromNode
    if((name %in% conduits_sf$ToNode) == F & length(which(conduits_sf$FromNode == name)) > 1){
      junctions$tag[junctions$Name %in% name] <- "hill" # "break"
    }

  }
  return(junctions)
}

#' correct artificial outfalls
#' @keywords internal
correct_artificial_outfalls <- function(junctions, conduits_sf){
  # select artifical outfalls:
  outfall_nodes_art <- junctions$Name[junctions$tag == "outfall_artificial"]

  for(node in outfall_nodes_art){
    # correct outfall node: transform outfall to start
    FromNode <- conduits_sf$FromNode[conduits_sf$ToNode == node]

    # adjust depths
    junctions$Bottom[junctions$Name == node] <- junctions$Bottom[junctions$Name == FromNode] + 0.01 # plus 1 mm

    # switch from and to node
    conduits_sf$FromNode[conduits_sf$ToNode == node] <- node
    conduits_sf$ToNode[conduits_sf$ToNode == node] <- FromNode

  }
  return(list(junctions = junctions, conduits = conduits_sf))
}

#' correct sinks
#' @keywords internal
correct_sinks <- function(junctions, conduits_sf, outfalls, short_cut_sinks, direct_drainage_sinks, ds){

  sink_at_end <- NULL
  count <- 1
  while(TRUE){

    # add tag "to_outfall"
    to_outfall_nodes <- conduits_sf$FromNode[conduits_sf$ToNode %in% outfalls$Name]
    junctions$tag[junctions$Name %in% to_outfall_nodes] <- "to_outfall"


    if(!("sink" %in% junctions$tag)){
      break
    }
    sink <- junctions$Name[junctions$tag == "sink"]

    for(sink in sink){
      # calculate depth of local sinks
      to_sink <- conduits_sf$FromNode[conduits_sf$ToNode == sink]

      to_type <- NULL
      for(n in to_sink){
        to_type <- c(to_type, junctions$tag[junctions$Name == n])
      }


      # drain towards "to_outfall"
      if("to_outfall" %in% to_type){
        to_outfall <- to_sink[to_type == "to_outfall"]
        name_conduit <- conduits_sf$Name[conduits_sf$FromNode == to_outfall & conduits_sf$ToNode == sink]
        conduits_sf$FromNode[conduits_sf$Name == name_conduit] <- sink
        conduits_sf$ToNode[conduits_sf$Name == name_conduit] <- to_outfall

      }else{

        # keep only nodes which are not "start" nodes in to_sink:
        to_sink <- to_sink[to_type != "start"]

        # keep nodes that were not sinks at the end the network (avoid creating artificial outfalls)
        to_sink <- to_sink[to_sink != sink_at_end]


        if(length(to_sink) > 1){
          neighbours <- junctions$Bottom[junctions$Name %in% to_sink]

          if(any(diff(neighbours) < ds, length(neighbours) == 1) & junctions$Top[junctions$Name %in% sink] > mean(neighbours)){
            # ds = depression storage (limit to either correct sinks or add a new conduit)

            # set bottom height of sink to mean height of neighbouring junctions
            junctions$Bottom[junctions$Name %in% sink] <- mean(neighbours)


            # define flow directions of conduits to sink:
            for(from in to_sink){
              name_conduit <- conduits_sf$Name[(conduits_sf$ToNode == sink) & (conduits_sf$FromNode == from)]
              if(junctions$Bottom[junctions$Name == from] < junctions$Bottom[junctions$Name == sink]){
                conduits_sf$ToNode[conduits_sf$Name == name_conduit] <- from
                conduits_sf$FromNode[conduits_sf$Name == name_conduit] <- sink
              }else{
                conduits_sf$ToNode[conduits_sf$Name == name_conduit] <- sink
                conduits_sf$FromNode[conduits_sf$Name == name_conduit] <- from
              }
            }

          }else{
            if(any(diff(neighbours) > ds, junctions$Top[junctions$Name %in% sink] < mean(neighbours))){

              # # add a new conduit to ensure the permeability of the network:
              if(short_cut_sinks){
                distances_sink <- as.numeric(sf::st_distance(junctions[junctions$Name == sink,], junctions))
                distances_sink <- data.frame(distances = distances_sink, Name = junctions$Name, Top = junctions$Top,
                                             stringsAsFactors = F)

                distances_sink <- distances_sink[distances_sink$distance > 0 & !distances_sink$Name %in% to_sink &
                                                   distances_sink$Top < junctions$Top[junctions$Name == sink],]

                name_nearest_lower_junction <-
                  distances_sink$Name[distances_sink$distances == min(distances_sink$distances)]


                new_conduit <- rbind(junctions[junctions$Name == name_nearest_lower_junction,], junctions[junctions$Name == sink,])
                new_conduit <- dplyr::summarize(new_conduit, .groups = "drop_last")
                new_conduit <- sf::st_cast(new_conduit, "LINESTRING")

                new_conduit$ToNode <- name_nearest_lower_junction

              }

              # add a connection to the nearest outfall:
              if(direct_drainage_sinks){

                name_nearest_outfall <- as.numeric(sf::st_distance(junctions[junctions$Name == sink,], outfalls))
                name_nearest_outfall <- outfalls$Name[which.min(name_nearest_outfall)]

                # swmm does not allow two links at one outfall; therefore, add another outfall:
                new_outfall_name <- paste0("n_",name_nearest_outfall)
                outfalls <- rbind(outfalls, outfalls[outfalls$Name == name_nearest_outfall,])
                outfalls[nrow(outfalls),"Name"] <- new_outfall_name
                outfalls[nrow(outfalls),"geometry"] <- outfalls[nrow(outfalls),"geometry"] + runif(1, min = 0, max = 2)

                new_conduit <- rbind(outfalls[nrow(outfalls),"Name"], junctions[junctions$Name == sink,"Name"])
                new_conduit <- dplyr::summarize(new_conduit, .groups = "drop_last")
                new_conduit <- sf::st_cast(new_conduit, "LINESTRING")
                new_conduit$ToNode <- new_outfall_name

                count <-  count +1
              }

              new_conduit$Name <- paste0("nC_",nrow(conduits_sf))
              new_conduit$FromNode <- sink
              new_conduit$L1 <- nrow(conduits_sf) + 1

              conduits_sf <- rbind(new_conduit, conduits_sf)

            }

          }
        }else{
          link <- conduits_sf$Name[conduits_sf$FromNode == to_sink & conduits_sf$ToNode == sink]
          conduits_sf$FromNode[conduits_sf$Name == link] <- sink
          conduits_sf$ToNode[conduits_sf$Name == link] <- to_sink

          # remember conduit and sink node
          sink_at_end <- sink

        }
      }



    }

    # check for sinks again
    junctions <- error_checks(junctions, conduits_sf, outfalls)

  }

  return(list(junctions = junctions, conduits_sf = conduits_sf, outfalls = outfalls))
}

#' track path from outfall to crossing
#' @keywords internal
linear_path <- function(next_to, outfalls, junctions, conduits_sf){

  path_slope_error <- list()

  while(TRUE){

    tag <- junctions$tag[junctions$Name == next_to]
    path_slope_error[[next_to]] <- tag

    if(tag %in% list("start", "crossing", "loop", "hill"))
    {
      break
    }

    next_to <- conduits_sf$FromNode[conduits_sf$ToNode == next_to]

  }

  return(list(next_to, tag, path_slope_error))
}

#' track path from junction (nearest to outfall) to start
#' @keywords internal
backtracking_paths <- function(start, outfalls, junctions, conduits_sf, next_to_seen = list()){

  triple <- linear_path(start, outfalls, junctions, conduits_sf)
  last_to <- triple[[1]]
  last_to_tag <- triple[[2]]
  path <- triple[[3]]

  paths <- list()
  skipped_links <- list()

  if(last_to_tag == "crossing"){

    for(next_to in conduits_sf$FromNode[conduits_sf$ToNode == last_to]){

      if(next_to %in% names(next_to_seen) & junctions$tag[junctions$Name == next_to] == "crossing"){

        if(next_to_seen[[next_to]] >= 10){
          skipped_links[[start]] <- c(start, next_to)

          next
        }

      }

      if(is.null(next_to_seen[[next_to]])){
        next_to_seen[[next_to]] <- 1
      }else{
        next_to_seen[[next_to]] <- next_to_seen[[next_to]] + 1
      }


      list_last_call <- backtracking_paths(next_to, outfalls, junctions, conduits_sf, next_to_seen)
      segments <- list_last_call$paths
      next_to_seen <- list_last_call$next_to_seen
      skipped_links <- c(list_last_call$skipped_links, skipped_links)

      for(segment in segments){
        extended_path <- c(path, segment)
        paths <- c(paths, list(extended_path))
      }

    }
  }else{
    paths <- list(path)
  }

  return(list(paths = paths, next_to_seen = next_to_seen, skipped_links = skipped_links))
}

#' track paths from all outfalls to start
#' @keywords internal
track_outfall_to_start_connectivity <- function(outfalls, junctions, conduits_sf){

  list_tracks <- list()
  for(i in 1:nrow(outfalls)){
    to_outfall <- conduits_sf$FromNode[conduits_sf$ToNode == outfalls$Name[i]]
    list_tracks[[to_outfall]] <- backtracking_paths(to_outfall, outfalls, junctions, conduits_sf)$paths
  }

  return(list_tracks)
}

#' calculate slopes
#' @keywords internal
calculateSlope <- function(junctions, conduits_sf, outfall){

  if(!any(colnames(conduits_sf) %in% "OutOffset")){
    warning("Column OutOffset is missing in conduits_sf")
  }
  if(!any(colnames(conduits_sf) %in% "InOffset")){
    warning("Column OutOffset is missing in conduits_sf")
  }

  # calculate slope
  conduits_sf$Slope <- NA

  for (i in 1:length(conduits_sf$Slope)){
    if(conduits_sf$ToNode[i] %in% outfall$Name){
      conduits_sf$Slope[i] <- (junctions$Bottom[junctions$Name == conduits_sf$FromNode[i]] - outfall$Elevation[outfall$Name == conduits_sf$ToNode[i]])/conduits_sf$Length[i]
    }else{
      conduits_sf$Slope[i] <- ((junctions$Bottom[junctions$Name == conduits_sf$FromNode[i]] + conduits_sf$InOffset[i]) - (junctions$Bottom[junctions$Name == conduits_sf$ToNode[i]] + conduits_sf$OutOffset[i]))/conduits_sf$Length[i]
    }
  }

  return(conduits_sf$Slope)
}

#' calculate new bottom height
#' @keywords internal
calculate_new_bottom_height <- function(length, slope, FromHeight){
  (((slope*length)) - FromHeight) * (-1)
}

#' slope and junction depth adjustment algorithm
#' @keywords internal
adjustSlopesAndJuncDepths <- function(conduits, junctions, outfall, min_slope, max_slope, min_junc_depth, mean_junc_depth, max_junc_depth){

  # prepare datasets to store adjustments:
  junctions$Bottom <- junctions$Top - mean_junc_depth
  junctions$max_bottom <- junctions$Top - max_junc_depth
  junctions$min_bottom <- junctions$Top - min_junc_depth

  # calculate slopes to be corrected:
  conduits$Slope <- calculateSlope(junctions, conduits, outfall)

  # store paths:
  network_paths <- track_outfall_to_start_connectivity(outfall, junctions, conduits)

  # return massage for cases when a stackoverflow occurs:
  message(" ... backtracking_paths(): if a stack overflow occurs, delete loops in the street network manually or set either break_closed_loops == TRUE or break_loops == TRUE")

  for(j in 1:length(network_paths)){
    for(i in 1:length(network_paths[[j]])){
      connections <- names(network_paths[[j]][[i]])
      junc <- length(connections)

      start <- dplyr::last(connections)

      while(TRUE){
        if((junc - 1) == 0){
          break
        }

        conduits$Slope <- calculateSlope(junctions, conduits, outfall)

        from <- connections[junc] #start
        to <- connections[junc-1] #end
        length <- conduits$Length[conduits$FromNode == from & conduits$ToNode == to]
        FromHeight <- junctions$Bottom[junctions$Name == from]
        slope <- conduits$Slope[conduits$FromNode == from & conduits$ToNode == to]
        minDepthTo <- junctions$min_bottom[junctions$Name == to]
        maxDepthTo <- junctions$max_bottom[junctions$Name == to]

        # a) normal case no changes needed: min_slope < slope & slope < max_slope){

        # b) if the slope is not between min and max slope
        if(slope < min_slope | slope > max_slope){

          #ToHeight <- junctions$Bottom[junctions$Name == to]

          if(slope < min_slope){
            ToHeight <- calculate_new_bottom_height(length, min_slope, FromHeight)

            if(ToHeight < maxDepthTo | ToHeight > minDepthTo){

              if(ToHeight < maxDepthTo){
                # bei flachem Gefälle zu tiefe junction
                message(paste(" ... from ", from,"- to", to, "slope < min_slope & ToHeight < maxDepthTo"))
              }

              if(ToHeight > minDepthTo){
                # enlarge slope to meet minDepthTo
                ToHeight <- minDepthTo

                newSlope <- (FromHeight - ToHeight)/length

                if(newSlope > max_slope){
                  message (paste(" ... from ", from,"- to", to,  "slope < min_slope & ToHeight > minDepthTo & newSlope > max_slope"))
                }

              }

            }



          }

          if(slope > max_slope){
            ToHeight <- calculate_new_bottom_height(length, max_slope, FromHeight)

            if(ToHeight < maxDepthTo | ToHeight > minDepthTo){
              if(ToHeight < maxDepthTo){
                message(paste("... from ", from,"- to", to,  "slope > max_slope & ToHeight < maxDepthTo"))
              }
              if(ToHeight > minDepthTo){
                ToHeight <- minDepthTo
                message(paste("... from ", from,"- to", to,  "slope > max_slope & ToHeight > minDepthTo"))
              }
            }
          }


          # take the smaller value of adjusted and original height
          junctions$Bottom[junctions$Name == to] <- min(c(ToHeight, junctions$Bottom[junctions$Name == to]))

          # remember the bottom height of the current pipe to calculate Offsets later:
          conduits$OutOffset[conduits$FromNode == from & conduits$ToNode == to] <- ToHeight

        }else{
          # no slope adjustments needed (slope between min and max)
          ToHeight <- junctions$Bottom[junctions$Name == to]

          # take the smaller value of adjusted and original height
          junctions$Bottom[junctions$Name == to] <- min(c(ToHeight, junctions$Bottom[junctions$Name == to]))

          # remember the bottom height of the current pipe to calculate Offsets later:
          conduits$OutOffset[conduits$FromNode == from & conduits$ToNode == to] <- ToHeight

        }

        junc <- junc - 1

      }

    }
  }

  # calculate final Offsets:
  for( i in 1:length(conduits$Name)){
    if(all(conduits$ToNode[i] != outfall$Name)){
      conduits$OutOffset[i] <- conduits$OutOffset[i] - junctions$Bottom[junctions$Name == conduits$ToNode[i]]
    }
  }

  # adjust slope outfall
  for(i in 1:nrow(outfall)){
    to <- outfall$Name[i]
    from <- conduits$FromNode[conduits$ToNode == to]

    length <- conduits$Length[conduits$FromNode == from & conduits$ToNode == to]
    FromHeight <- junctions$Bottom[junctions$Name == from]
    ToHeight <- calculate_new_bottom_height(length, min_slope, FromHeight)
    outfall$Elevation[i] <- ToHeight
  }

  # plot final cross sections:
  for(j in 1:length(network_paths)){
    for(i in 1:length(network_paths[[j]])){
      connections <- names(network_paths[[j]][[i]])
      if(length(connections) <= 1){
        next
      }
      plot_cross_section_junctions(junctions, conduits, connections)
    }
  }

  return(list(junctions = junctions, conduits = conduits, outfalls = outfall))

}

#' breaks in network
#' @keywords internal
breaks_in_network <- function(break_closed_loops, break_loops, breaks_at_hills, delete_disconnected, junctions, outfalls, conduits_sf){
  if(break_closed_loops){
    #### if an stack overflow occurs because of closed loops in the network devide network as follows: ####
    # remember skipped conduits:
    skipped_links_all <- data.frame(start = NULL, skipped = NULL, stringsAsFactors = F)

    for(i in 1:nrow(outfalls)){
      to_outfall <- conduits_sf$FromNode[conduits_sf$ToNode == outfalls$Name[i]]
      skipped_links_all_new <- backtracking_paths(to_outfall, outfalls, junctions, conduits_sf)$skipped_links
      skipped_links_all_new <- do.call(rbind, skipped_links_all_new)
      skipped_links_all <- rbind(skipped_links_all, skipped_links_all_new)
    }

    if(length(skipped_links_all$start) > 0){
      colnames(skipped_links_all) <- c("start", "skipped")
      skipped_links_all$start <- as.character(skipped_links_all$start)
      skipped_links_all$skipped <- as.character(skipped_links_all$skipped)

      # find link between start and skipped
      delete_links <- c()
      for(i in 1:nrow(skipped_links_all)){
        if(any(conduits_sf$FromNode == skipped_links_all$skipped[i] & conduits_sf$ToNode == skipped_links_all$start[i])){
          delete_links_new  <- which(conduits_sf$FromNode == skipped_links_all$skipped[i] & conduits_sf$ToNode == skipped_links_all$start[i])
          delete_links <- c(delete_links, delete_links_new)
        }else{
          skipped <- linear_path(skipped_links_all$start[i], outfalls, junctions, conduits_sf)[[1]]
          delete_links_new <- which(conduits_sf$FromNode == skipped & conduits_sf$ToNode == skipped_links_all$start[i])
          delete_links <-  c(delete_links, delete_links_new)
        }
      }

      #  exclude links for slope correction
      conduits_sf <- conduits_sf[-delete_links,]

      # new tags:
      junctions <- error_checks(junctions, conduits_sf, outfalls)

      # plot tags...
      plot_tagged_junctions(junctions, outfalls, conduits_sf, main = "tags after deletion of conduits")

      # repeat sink correction:
      list_sink_corrections <- correct_sinks(junctions, conduits_sf, outfalls)
      junctions <- list_sink_corrections$junctions
      conduits_sf <- list_sink_corrections$conduits_sf

      # ... and plot:
      plot_tagged_junctions(junctions, outfalls, conduits_sf, main = "repeated correction of sinks")


    }
  }

  if(breaks_at_hills){

    conduits_sf$Length <- as.numeric(sf::st_length(conduits_sf))

    while("hill" %in% junctions$tag){
      #### delete shorter conduit at hill:
      for(j in junctions$Name[junctions$tag == "hill"]){
        conduits_from_hill <- conduits_sf$Name[conduits_sf$FromNode %in% j]
        index_skip <- which.min(conduits_sf$Length[conduits_sf$Name %in% conduits_from_hill])
        conduits_sf <- conduits_sf[!conduits_sf$Name==conduits_from_hill[index_skip],]
      }

      junctions <- error_checks(junctions, conduits_sf, outfalls)
    }

    # delete outfalls that are disconnected now:
    for(out in outfalls$Name){
      if(out %in% conduits_sf$ToNode){
        next
      }else{
        outfalls <- outfalls[outfalls$Name != out,]
      }
    }


    # plot tags...
    plot_tagged_junctions(junctions, outfalls, conduits_sf, main = "breaks at hills")
    message("a. network is devided at hill nodes")


  }

  if(delete_disconnected){
    # plot conduits that are not connected to an outfall:
    junctions <- error_checks(junctions, conduits_sf, outfalls)
    plot_disconnected_conduits(junctions, conduits_sf, outfalls, col_arrow = "red")

    list_disconnected <- delete_disconnected_network_parts(junctions, conduits_sf, outfalls)
    junctions <- list_disconnected$junctions
    conduits_sf <- list_disconnected$conduits

    j_coords <- as.data.frame(sf::st_coordinates(junctions))
    j_coords$Name <- junctions$Name
    plot_drainage_network(j_coords, outfalls, conduits_sf, main = "final network after deletion of disconnected parts", col_arrow = "black")
  }

  if(break_loops){
    #### delete conduit at junctions which drain in two directions:
    junctions <- error_checks(junctions, conduits_sf, outfalls)

    for(name in unique(junctions$Name)){
      if(all(length(which(conduits_sf$FromNode == name)) > 1 & length(which(conduits_sf$ToNode == name)) >= 1 &
             !conduits_sf$ToNode[conduits_sf$FromNode == name] %in% outfalls$Name)){
        junctions$tag[junctions$Name == name] <- "loop"
      }
    }

    plot_tagged_junctions(junctions, outfalls, conduits_sf, main = "highlight loop initiating crossing")

    while(TRUE){

      if(!"loop" %in% junctions$tag){
        break
      }

      # track all paths
      all_paths <- track_outfall_to_start_connectivity(outfalls, junctions, conduits_sf)
      loop_nodes <- junctions$Name[junctions$tag == "loop"]

      # find paths to loop nodes and calculate path length:
      list_loop_paths <- vector("list", length(loop_nodes))
      names(list_loop_paths) <- loop_nodes

      for(node in loop_nodes){

        list_loop_paths[[node]] <-  list()
        list_loop_paths[[node]][["m"]] <-  NA
        list_loop_paths[[node]][["c"]] <- NA

        i = 1

        for(catchment in 1:length(all_paths)){
          for(path in 1:length(all_paths[[catchment]])){

            # proceed only if the loop node occurs on the path
            if(!node %in% names(all_paths[[catchment]][[path]])){
              next
            }

            nodes <- length(all_paths[[catchment]][[path]])
            node_names <- names(all_paths[[catchment]][[path]])

            # calculate path length:
            meter <- 0
            for(j in 2:nodes){
              L <- conduits_sf$Length[conduits_sf$ToNode == node_names[j-1] & conduits_sf$FromNode == node_names[j]]
              meter <- sum(meter, L)
            }

            # connecting conduit name:
            conduit_to_loop <- conduits_sf$Name[conduits_sf$ToNode == node_names[length(node_names)-1] & conduits_sf$FromNode == node]

            if(length(conduit_to_loop)>1){
              break
            }

            list_loop_paths[[node]][["m"]][i] <- meter
            list_loop_paths[[node]][["c"]][i] <- conduit_to_loop

            i = i +1

          }

        }


      }

      # keep shortest path to outfall:
      for(node in names(list_loop_paths)){
        index_m_min <- which.min(list_loop_paths[[node]]$m)
        c_del <- list_loop_paths[[node]]$c[-index_m_min]

        # delete conduits that are not the shortest path to the outfall:
        conduits_sf <- conduits_sf[!conduits_sf$Name %in% c_del,]
      }

      # tag again:
      junctions <- error_checks(junctions, conduits_sf, outfalls)

      for(name in unique(junctions$Name)){
        if(all(length(which(conduits_sf$FromNode == name)) > 1 & length(which(conduits_sf$ToNode == name)) >= 1 &
               !conduits_sf$ToNode[conduits_sf$FromNode == name] %in% outfalls$Name)){
          junctions$tag[junctions$Name == name] <- "loop"
        }
      }

    }

    plot_tagged_junctions(junctions, outfalls, conduits_sf, main = "network without loops")

    message("b. network is devided at loop initiating nodes")

  }

  return(list(conduits = conduits_sf, junctions = junctions, outfalls = outfalls))
}

#' delete junctions and conduits that are not connected to an outfall
#' @keywords internal
delete_disconnected_network_parts <- function(junctions, conduits, outfalls){

  # exclude outfall with only one draining junction connected:
  outfalls_con <- outfalls
  for(out in outfalls$Name){
    from <- conduits$FromNode[conduits$ToNode == out]
    if(any(conduits$ToNode == from | length(which(conduits$FromNode == from)) > 1)){
      next
    }else{
      outfalls_con <- outfalls_con[outfalls_con$Name != out,]

    }
  }

  # paths to outfall
  network_paths <- track_outfall_to_start_connectivity(outfalls_con, junctions, conduits)

  all_conduits_connected <- NULL
  for(j in 1:length(network_paths)){
    for(i in 1:length(network_paths[[j]])){
      names_junc <- names(network_paths[[j]][[i]])

      conduits_connected <- NULL
      for(n in 2:length(names_junc)){
        to <- names_junc[n-1]
        from <- names_junc[n]
        conduits_connected[n-1] <- conduits$Name[conduits$FromNode == from & conduits$ToNode == to]
      }


      if(is.null(all_conduits_connected)){
        all_conduits_connected <- conduits_connected
      }else{
        all_conduits_connected <- c(all_conduits_connected, conduits_connected)
      }
    }
  }

  conduits_disconnected <- conduits$Name[!(conduits$Name %in% unique(all_conduits_connected) | conduits$ToNode %in% outfalls$Name)]

  # filter for pipes that are connected:
  conduits <- conduits[!conduits$Name %in% conduits_disconnected,]

  # filter for junctions that are connected:
  junctions <- junctions[junctions$Name %in% conduits$FromNode | junctions$Name %in% conduits$ToNode,]


  return(list(junctions = junctions, conduits = conduits))


}

#' Create Thiessen polygons
#'
#' Simplified approach to define recharging surfaces for SWMM model.
#'
#' @param points Centroid points for Thiessen tesselation. Here junction nodes.
#' @param boundary_polygon model boundary polygon to cut Thiessen polygons.
#' @return A sf object containing Thiessen polygons.
#' @export
#' @rdname thiessenpolygons
thiessenpolygons <- function(points, boundary_polygon){

  # centroids
  coords <- sf::st_coordinates(points)

  # calculate thiessen polygons
  analyses <- deldir::deldir(coords[,"X"], coords[,"Y"])
  thiessen_list <- deldir::tile.list(analyses)

  # transform to sf
  polygons <- NA

  for( i in seq(along = thiessen_list)){

    polygon_coords <- cbind(thiessen_list[[i]]$x, thiessen_list[[i]]$y)
    polygon <- list(rbind(polygon_coords, polygon_coords[1,]))
    polygon <- sf::st_geometry(sf::st_polygon(polygon))

    if(any(is.na(polygons))){
      polygons <- polygon
    }else{
      polygons <- c(polygons, polygon)
    }

  }

  thiessen_polygons <- sf::st_sf(data.frame(ID = seq(along = thiessen_list), geometry = polygons))

  # add coordinate reference system
  sf::st_crs(thiessen_polygons) <- sf::st_crs(points)

  # reference to point features
  thiessen_polygons$centroid <- points$Name

  # cut with boundary polygon if avaiable and handle multipolygons
  if(!is.null(boundary_polygon)){
    thiessen_polygons <- suppressWarnings(sf::st_intersection(thiessen_polygons, boundary_polygon))
    thiessen_polygons <- sf::st_cast(thiessen_polygons, "MULTIPOLYGON")
    thiessen_polygons <- sf::st_cast(thiessen_polygons, "POLYGON", warn = F)

    for(c in unique(thiessen_polygons$centroid)){
      sub <- which(thiessen_polygons$centroid == c)
      if(length(sub) > 1){
        keep <- sf::st_distance(points[points$Name == c,], thiessen_polygons[sub,])
        keep <- which.min(as.numeric(keep))

        thiessen_polygons <- thiessen_polygons[-sub[-keep],]

      }
    }

  }

  # calculate areas
  thiessen_polygons$area <- as.numeric(sf::st_area(thiessen_polygons))

  # plot result:
  graphics::plot(sf::st_geometry(thiessen_polygons), border = "black", col = "white", axes = T, main = "Thiessen polygons with centroids")
  graphics::plot(sf::st_geometry(points), add = T, col = "red")

  return(thiessen_polygons[,c("ID","centroid", "area")])

}

#' combine thiessen polygons and land use information
#' @keywords internal
#' @param thiessen_polygons A sf object containing Thiessen polygons with junctions as centroids. Run function \link[urbandrain]{thiessenpolygons}.
#' @param landuse_sf A sf object containing polygons with land use informations as attributes.
#' @param landuse_classes A data frame specifying the percentage of imperviousness for landuse_sf.
#' @param re Alternative to landuse_sf: A raster file containing landuse information. For example derived from RapidEye Satelite Imagery.
#' @param area_re_pixel A data frame adding attributes to the re raster pixel values.
#' @param re_classes A data frame specifying the percentage of imperviousness for area_re_pixel.
#' @param width not implemented yet
#' @param slope not implemented yet
#' @param soil not implemented yet
#' @return A sf object including all information needed for the subcatchment input section in SWMM's input file.
#' @export
#' @rdname thiessen_to_subcatchments
thiessen_to_subcatchments <- function(landuse_sf = NULL, landuse_classes = NULL, re = NULL, area_re_pixel = NULL, re_classes = NULL, thiessen_polygons, width = NULL, slope = NULL, soil = NULL){

  if(!is.null(re)){
    if(all(!is.null(area_re_pixel) & !is.null(re_classes))){
      # if landcover data is given as RapidEye Raster data, calculate percent imperviousness as follows:
      # zonal statistics
      clip <- raster::extract(re, thiessen_polygons)
      clip <- lapply(clip, table)

      # impervious re classes
      impervious <- re_classes[-which(re_classes$perviousness == "pervious"),]

      thiessen_polygons$PercImperv <- NA

      for(id in 1:length(clip)){
        thiessen_polygons$PercImperv[id] <- as.numeric((sum(clip[[id]][names(clip[[id]]) %in% impervious$class]) * 100)/
                                                         sum(clip[[id]]))
      }
    }else{
      warning("Please specify the arguments area_re_pixel and re_classes")
    }

  }

  if(!is.null(landuse_sf)){
    if(!is.null(landuse_classes)){

      # intersect data:
      landuse_thiessen <- sf::st_intersection(landuse_sf, thiessen_polygons)
      # calculate polygon area of intersected data:
      landuse_thiessen$area_landuse <- as.numeric(sf::st_area(landuse_thiessen))

      # add a factor for impervious surfaces as defined in landuse_classes:
      landuse_thiessen$factor_imperv <- NA
      for(i in landuse_classes$class){
        landuse_thiessen$factor_imperv[landuse_thiessen$landuse == i] <- landuse_classes$imperviousness[landuse_classes$class == i]
      }

      # calculate percentage of impervious area per thiessen polygon:
      for(id in thiessen_polygons$ID){
        sel <- landuse_thiessen[landuse_thiessen$ID == id,]
        thiessen_polygons$PercImperv[thiessen_polygons$ID == id] <- sum(sel$area_landuse * sel$factor_imperv)/sel$area[1] * 100
      }

    }else{
      warning("Please specify the argument landuse_classes")
    }

  }

  # new subcatchment name
  thiessen_polygons$Name <- paste0("S_",thiessen_polygons$ID)

  # Area in ha
  thiessen_polygons$Area <- thiessen_polygons$area/10000

  # RouteTo: pervious and impervious subareas route to outlet node
  thiessen_polygons$RouteTo <- "OUTLET"
  thiessen_polygons$Outlet <- thiessen_polygons$centroid

  # (PctRouted)--> optional
  thiessen_polygons$PctRouted <- 100

  if(is.null(width)){
    # ideal rectangular --> A_T = 2*W
    thiessen_polygons$Width <- thiessen_polygons$area/2
  }else{
    message("flow path analysis is not implemented yet...")
  }

  if(is.null(slope)){
    # default slope value
    thiessen_polygons$Slope <- 1
  }else{
    # calculate mean slope from dtm
  }

  if(is.null(soil)){
    # default soil value
    thiessen_polygons$Soil <- "Sand"
  }else{
    # add soil information?
    message("intersection with soil maps is not implemented yet ...")
  }

  subcatchment <- thiessen_polygons[,c("Name", "Outlet", "Area", "RouteTo", "Soil", "PercImperv", "Width", "Slope", "PctRouted")]

  return(subcatchment)
}

#' Pipe diameter sizing algorithm
#'
#' The conduit diameters are calculated and designed for a user defined storm event. Initially the sizing algorithm runs a SWMM model with all conduits having a diameter of 0.3 m by default (or smaller if defined). After each iteration, it is checked whether conduits still surcharge. If so, diameters of the surcharging conduits are increased by a nominal diameter. The algorithm stops if no conduit surcharge occurs any more. The final inp file is stored.
#'
#' @param inp A object of class inp. Containing a complete SWMM model.
#' @param path_out A directory path to write the results to.
#' @param conduits_sf A sf object containing polyline information.
#' @param junctions A sf object containing junction point information.
#' @param outfalls A sf object containing outfall point informations.
#' @param target Either node_flooding or conduit_surcharge.
#' @param DN_start Initial pipe diameter (unit: mm).
#' @param rpt_filename File name for SWMM's report file.
#' @param inp_filename File name for SWMM's input file.
#' @param link_filename File name for a shp file containing the sized conduits.
#' @return A sf object containing polylines with all attributes needed for the conduit section in SWMM's input file.
#' @export
#' @rdname pipeSizingAlgorithm
pipeSizingAlgorithm <- function(inp, path_out, conduits_sf, junctions, outfalls, target, DN_start = 300, rpt_filename = "artificial_SWMM_model.rpt", inp_filename = "artificial_SWMM_model.inp", link_filename = "links_artificial_SWMM_format.shp", swmm5_exec){

  # Nennweiten [m]:
  DN_mm <- c(100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1800, 2000,
             2200, 2500, 2800, 3000, 3200, 3500, 4000)  #100,150, 200, 250,
  DN_m <- DN_mm/1000

  DN_start_index <- which(DN_mm == DN_start)
  steps <- NULL
  for (i in DN_start_index:(length(DN_m)-1)){
    step <- DN_m[i+1] - DN_m[i]
    steps <- c(steps, step)
  }

  if(target == "conduit_surcharge"){
    s <- 1
    while(TRUE){
      # initiate model runs as long as surchage occurs
      rpt <- swmmr::read_rpt(file.path(path_out, rpt_filename))

      if(is.null(rpt$conduit_surcharge_summary$Conduit)){
        break
      }

      # surcharging conduits:
      conduit_full <- rpt$conduit_surcharge_summary$Conduit

      # enlarge selected pipes:
      inp$xsections$Geom1[inp$xsections$Link %in% conduit_full] <- inp$xsections$Geom1[inp$xsections$Link %in% conduit_full] + steps[s]

      # plot new diameters:
      conduits_sf$Geom1 <- inp$xsections$Geom1
      plot_drainage_network_pipe_diameters(junctions, outfalls, conduits_sf, main = "pipe diameters during pipe sizing", DN_m)

      # save and run the model
      swmmr::write_inp(inp, file.path(path_out, inp_filename))
      swmmr::run_swmm(file.path(path_out, inp_filename),exec = swmm5_exec)

      # next pipe diameter
      s <- s + 1

    }
  }

  if(target == "node_flooding"){
    s <- 1
    while(TRUE){
      # initiate model runs as long as surchage occurs
      rpt <- swmmr::read_rpt(file.path(path_out, rpt_filename))

      if(is.null(rpt$node_flooding_summary$Node)){
        break
      }

      # flooded nodes:
      node_flooded <- rpt$node_flooding_summary$Node

      # conduits to enlarge:
      conduit_full <- NULL
      for(nf in node_flooded){
        c1 <- inp$conduits$Name[inp$conduits$FromNode == nf]
        c2 <- inp$conduits$Name[inp$conduits$ToNode == nf]
        conduit_full <- c(conduit_full, c1, c2)
      }
      conduit_full <- unique(conduit_full)


      # enlarge selected pipes:
      inp$xsections$Geom1[inp$xsections$Link %in% conduit_full] <- inp$xsections$Geom1[inp$xsections$Link %in% conduit_full] + steps[s]

      # plot new diameters:
      conduits_sf$Geom1 <- inp$xsections$Geom1
      plot_drainage_network_pipe_diameters(junctions, outfalls, conduits_sf, main = "pipe diameters during pipe sizing", DN_m)

      # save and run the model
      swmmr::write_inp(inp, file.path(path_out, inp_filename))
      swmmr::run_swmm(file.path(path_out, inp_filename),exec = swmm5_exec)

      # next pipe diameter
      s <- s + 1


    }
  }


  # adjust pipe diameter in line.shp
  line_file <- file.path(path_out, link_filename)

  if(file.exists(line_file)){
    file.remove(file.path(path_out, list.files(path = path_out, pattern = gsub(".shp", "", link_filename))))
  }

  conduits_sf$Geom1 <- inp$xsections$Geom1
  sf::st_write(conduits_sf, file.path(path_out, link_filename))

  # plot new pipe diameters:
  plot_drainage_network_pipe_diameters(junctions, outfalls, conduits_sf, main = "pipe diameters after pipe sizing", DN_m)

  return(conduits_sf)
}
