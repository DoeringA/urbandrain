#' Create a drainage network including junctions, conduits and outfalls.
#'
#' Fully automated generation of a drainage network inclduding junctions, conduits and outfalls based on open source spatial data. Requires street polylines, a digital terrain model, outfall point data and some user defined parameters as input. Returns a list of the sf objects junctions, conduits and outalls that can be converted to a SWMM *.inp file using swmmr functions afterwards.
#'
#' @param streets Street polyline data of class sf. On basis of the street polylines the drainage network is setup. Street polyline data can be extracted for example from openstreetmaps.
#' @param dtm Digital terrain information either of class RasterLayer or as a triangulated grid of class sf. Surface height data is required to add heights to junctions and to define the flow directions in the conduits subsequently.
#' @param outfalls Outfall point data of class sf. The network drains towards these points.
#' @param crs_default A proj4string (if dtm is of class RasterLayer) or EPSG Code to specify the default coordinate reference system. The measure of the default crs must be in m.
#' @param buffer Minimum distance between street nodes (unit: meter). If distance is smaller street nodes are summarized. A parameter for quality checks on street polyline data.
#' @param snap_dist Minimum distance between a polyline and an ending node (unit: meter). If distance is smaller street polylines are connected. A parameter for quality checks on street polyline data.
#' @param epsilon Distance to ideal line between start and end point (unit: meter). Parameter from the Douglas Peucker Algorithm that is needed to simplify the shape of polylines. See \link[kmlShape]{DouglasPeuckerEpsilon} for further information.
#' @param lim The maximum distance between two junctions (unit: meter). Can be interpreted as maximum pipe length. To avoid very short pipes this values can be exceeded occasionally.
#' @param min_junc_depth Minimum junction depth (unit: meter).
#' @param mean_junc_depth Mean junction depth (unit: meter).
#' @param max_junc_depth Maximum junction depth (unit: meter).
#' @param min_slope Minimum slope value (unit: -)
#' @param max_slope Maximum slope value (unit: -)
#' @param ds Threshold for local sinks in the network (unit: meter). If the depth of a local sink is greater than ds (depression storage). A shortcut is implemented either to the nearest outfall if direct_drainage_sinks is set TRUE or to the nearest lower junction if short_cut_sinks is set TRUE.
#' @param stepwise Stepwise construction of the network if set to TRUE. When a stepwise construction is chosen, at first a base network including junctions at start, end and crossing fix points is computed. Then the local sink nodes are corrected.
#' @param break_closed_loops -
#' @param delete_disconnected -
#' @param break_loops If set to TRUE: The drainage network is forked at junctions that have no incoming but more than one outgoing pipe. Parameter to increase the linearity of the network.
#' @param breaks_at_hills If set to TRUE: The drainage network is forked at junctions that have incoming and more than one outgoing pipe. Parameter to increase the linearity of the network.
#' @param short_cut_sinks If set to TRUE: drainage of local sinks that are deeper than ds towards the nearest lower junction.
#' @param direct_drainage_sinks If set to TRUE (default): drainage of local sinks that are deeper than ds towards the nearest outfall node.
#' @return A list of sf objects junctions, conduits and outfalls.
#' @export
#' @rdname create_drainage_network
#' @importFrom stats runif

create_drainage_network <- function(streets, dtm, outfalls, crs_default , buffer, snap_dist, epsilon, lim,  min_junc_depth, mean_junc_depth, max_junc_depth,  min_slope, max_slope, ds, stepwise, break_closed_loops,  delete_disconnected, breaks_at_hills, break_loops,  short_cut_sinks, direct_drainage_sinks){

  # ... check input data for structure
  if(all(class(streets) != "sf")){
    stop("class(streets) != sf")
  }

  if(all(class(dtm) != "RasterLayer" & class(dtm) != "sf")){
    stop("class(dtm) != RasterLayer & class(dtm) != sf")
  }

  if(all(class(outfalls) != "sf")){
    stop("class(outfalls) != sf")
  }

  if(!any(is.numeric(c(buffer, snap_dist, epsilon, lim,  min_junc_depth, mean_junc_depth, max_junc_depth,  min_slope, max_slope, ds)))){
    stop("buffer, snap_dist, epsilon, lim,  min_junc_depth, mean_junc_depth, max_junc_depth,  min_slope, max_slope, ds have to be numeric")
  }

  if(!any(is.logical(c(stepwise, break_closed_loops,  delete_disconnected, breaks_at_hills, break_loops,  short_cut_sinks, direct_drainage_sinks)))){
    stop("stepwise, break_closed_loops,  delete_disconnected, breaks_at_hills, break_loops,  short_cut_sinks, direct_drainage_sinks only TRUE or FALSE are allowed")
  }

  # ... transform input data to crs_default
  streets <- sf::st_transform(streets, crs_default)

  outfalls <- sf::st_transform(outfalls, crs_default)

  if(any(class(dtm) == "sf")){
    dtm <- sf::st_transform(dtm, crs_default)
  }

  if(any(class(dtm) == "RasterLayer")){
    dtm <- suppressWarnings(raster::projectRaster(dtm, crs = crs_default))
  }


  if(!stepwise){
    #### 1. Create junctions along the street polylines: ####
    junctions <- create_junctions(streets = streets, buffer = buffer, snap_dist = snap_dist, epsilon = epsilon, lim = lim, junc_depth = NULL, dtm = dtm, crs_default = crs_default, pre_def_junctions = NULL, pre_def_conduits = NULL)
    message("1. junctions created")

    #### 2. Add junctions heights: ####
    junctions$Top <- add_junction_heights(junctions, dtm)
    # calculate bottom height
    junctions$Bottom <- junctions$Top - mean_junc_depth

    message("2. junction heights added to initial junctions")

    #### 3. Define conduits along the topography: ####
    network_list <- create_conduits_along_topography(vertices = junctions, main = "flow directions following the topography")
    junctions <- network_list$junctions
    conduits_sf <- network_list$conduits_sf

    # progress...
    message("3. flow directions in conduits are defined following the topography")

    #### 4. Add outfall conduit connection: ####
    network_list <- add_outfall_connection(outfalls = outfalls, junctions = junctions, conduits_sf = conduits_sf, min_slope = min_slope, junc_depth = min_junc_depth, lim = lim)
    junctions <- network_list$junctions
    conduits_sf <- network_list$conduits
    outfalls <- network_list$outfalls

    message("4. connections to outfalls are added")

    #### Tag junctions: ####
    junctions <- error_checks(junctions, conduits_sf, outfalls)

    # plot tags...
    plot_tagged_junctions(junctions, outfalls, conduits_sf, main = "add tags to junctions")

    #### 5. Correct artificial outfalls: ####
    corrected_outfalls <- correct_artificial_outfalls(junctions = junctions, conduits_sf = conduits_sf)
    junctions <- corrected_outfalls$junctions
    conduits_sf <- corrected_outfalls$conduits

    # add new tags:
    junctions <- error_checks(junctions, conduits_sf, outfalls)

    # ... and plot:
    plot_tagged_junctions(junctions, outfalls, conduits_sf, main = "corrected artificial outfalls")
    message("5. artificial outfalls are corrected and transformed to start nodes")

    #### Add name to conduits: ####
    conduits_sf$Name <- paste0("C_", 1:nrow(conduits_sf))

    #### 6. Correct sinks: ####
    list_sink_corrections <- correct_sinks(junctions, conduits_sf, outfalls, short_cut_sinks, direct_drainage_sinks, ds)
    junctions <- list_sink_corrections$junctions
    conduits_sf <- list_sink_corrections$conduits_sf
    outfalls <- list_sink_corrections$outfalls

    # ... and plot:
    junctions <- error_checks(junctions, conduits_sf, outfalls)
    plot_tagged_junctions(junctions, outfalls, conduits_sf, main = "corrected sinks")
    plot_short_cuts(junctions, conduits_sf, outfalls, "red")
    message("6. local sinks in network are corrected")

    #### 7. Devide network (optional): ####
    message("7. optional breaks in network to increase linearity:")
    if(any(break_closed_loops, break_loops, breaks_at_hills, delete_disconnected)){
      list_breaks <- breaks_in_network(break_closed_loops, break_loops, breaks_at_hills, delete_disconnected, junctions = junctions, outfalls = outfalls, conduits_sf = conduits_sf)
      junctions <- list_breaks$junctions
      conduits_sf <- list_breaks$conduits
      outfalls <- list_breaks$outfalls
    }

    #### 8. Check consistency of network from outfall to start or hill ####
    junctions <- error_checks(junctions, conduits_sf, outfalls)
    # plot paths:
    for (i in 1:nrow(outfalls)){
      plot_path_to_outfall(junctions, outfalls, oid = i, conduits_sf, main = paste0("path to outfall ", i), col_arrow = "red")
    }

    message("8. outfall catchments are plotted")

    #### 9. Slopes and junction depths are adjusted within the given range ####
    conduits_sf$Length <- sf::st_length(conduits_sf)
      conduits_sf$Length <-  as.numeric(conduits_sf$Length)
    conduits_sf$InOffset <- 0
    conduits_sf$OutOffset <- 0

    # calculate initial slopes and plot them...
    conduits_sf$Slope <- calculateSlope(junctions = junctions, conduits_sf = conduits_sf, outfall = outfalls)
    graphics::boxplot(conduits_sf$Slope, main = "conduit slopes before adjustment")
    graphics::boxplot(junctions$Top - junctions$Bottom, main = "junction depths before slope adjustment")
    plot_slopes(junctions = junctions, conduits_sf = conduits_sf, outfalls = outfalls, min_slope = min_slope, max_slope = max_slope, main = "slopes before adjustment")

    # correct slopes and junction depths in the range min/max slope and junction depth
    list_slope_adjustments <- adjustSlopesAndJuncDepths(conduits = conduits_sf, junctions = junctions, outfall = outfalls,
                                           min_slope = min_slope, max_slope = max_slope,
                                           min_junc_depth = min_junc_depth, max_junc_depth = max_junc_depth,
                                           mean_junc_depth = mean_junc_depth)

    outfalls <- list_slope_adjustments$outfall
    junctions <- list_slope_adjustments$junctions
    conduits_sf <- list_slope_adjustments$conduits

    # calculate final slopes ..
    conduits_sf$Slope <- calculateSlope(junctions = junctions, conduits_sf = conduits_sf, outfall = outfalls)
    # plot slopes ...
    graphics::boxplot(conduits_sf$Slope, main = "conduit slopes after adjustment")
    graphics::boxplot(junctions$Top - junctions$Bottom, main = "junction depths after slope adjustment")
    plot_slopes(junctions = junctions, conduits_sf = conduits_sf, outfalls = outfalls, min_slope = min_slope, max_slope = max_slope,main = "slopes after adjustment")

    message("9. Slopes and junction depths are adjusted within the given range (min_junc_depth, max_junc_depth, min_slope, max_slope)")

    #### 10. add remaining default parameters: ####
    conduits_sf$Shape <- "CIRCULAR"
    conduits_sf$Geom1 <- 0.3
    conduits_sf$Roughness <- 0.018

    message("10. missing default parameters are added")

    #### return ####
    return(list(junctions = junctions, conduits = conduits_sf, outfalls = outfalls))

  }else{
    #### 1. First call for creation of junctions (min possible amount of junctions): ####
    message("0. stepwise construction of the drainage network")

    junctions <- create_junctions(streets = streets, buffer = buffer, snap_dist = snap_dist, epsilon = 10000, lim = 10000, junc_depth = NULL, dtm = dtm, crs_default = crs_default, pre_def_junctions = NULL, pre_def_conduits = NULL)
    message("1. initial junctions are created")

    #### 2. Add junctions heights: ####
    junctions$Top <- add_junction_heights(junctions, dtm)
    # calculate bottom height
    junctions$Bottom <- junctions$Top - mean_junc_depth

    message("2. junction heights added to initial junctions")

    #### 3. Define conduits along the topography: ####
    network_list <- create_conduits_along_topography(vertices = junctions, main = "flow directions following the topography")
    junctions <- network_list$junctions
    conduits_sf <- network_list$conduits_sf

    # progress...
    message("3. flow directions in conduits are defined following the topography")

    #### 4. Add outfall conduit connection: ####
    network_list <- add_outfall_connection(outfalls = outfalls, junctions = junctions, conduits_sf = conduits_sf, min_slope = min_slope, junc_depth = min_junc_depth, lim = lim)
    junctions <- network_list$junctions
    conduits_sf <- network_list$conduits
    outfalls <- network_list$outfalls

    message("4. connections to outfalls are added")

    #### Tag junctions: ####
    junctions <- error_checks(junctions, conduits_sf, outfalls)

    # plot tags...
    plot_tagged_junctions(junctions, outfalls, conduits_sf, main = "add tags to junctions")

    #### 5. Correct artificial outfalls: ####
    corrected_outfalls <- correct_artificial_outfalls(junctions = junctions, conduits_sf = conduits_sf)
    junctions <- corrected_outfalls$junctions
    conduits_sf <- corrected_outfalls$conduits

    # add new tags:
    junctions <- error_checks(junctions, conduits_sf, outfalls)

    # ... and plot:
    plot_tagged_junctions(junctions, outfalls, conduits_sf, main = "corrected artificial outfalls")
    message("5. artificial outfall are corrected and tranformed to start nodes")

    #### Add name to conduits: ####
    conduits_sf$Name <- paste0("C_", 1:nrow(conduits_sf))

    #### 6. Correct sinks: ####
    list_sink_corrections <- correct_sinks(junctions, conduits_sf, outfalls, short_cut_sinks, direct_drainage_sinks, ds)
    junctions <- list_sink_corrections$junctions
    conduits_sf <- list_sink_corrections$conduits_sf

    # ... and plot:
    junctions <- error_checks(junctions, conduits_sf, outfalls)
    plot_tagged_junctions(junctions, outfalls, conduits_sf, main = "corrected sinks")
    plot_short_cuts(junctions, conduits_sf, outfalls, "red")

    message("6. local sinks in initial network were corrected")

    #### 7. Create junctions with user defined epsilon and lim values: ####
    network_list_2 <- create_junctions(streets, buffer, snap_dist, epsilon, lim, junc_depth = mean_junc_depth, dtm = dtm, crs_default, pre_def_junctions = junctions, pre_def_conduits = conduits_sf)

    junctions_u <- network_list_2$junctions
    conduits_u <- network_list_2$conduits

    # plot drainage network:
    coords_junctions <- sf::st_coordinates(junctions_u)
    coords_junctions <- as.data.frame(coords_junctions)
    coords_junctions$Name <- junctions_u$Name
    plot_drainage_network(coords = coords_junctions, outfalls = outfalls, conduits_sf= conduits_u, main = "user defined drainage network", col_arrow = "red")

    # Top height from dtm, bottom height interpolated from first junction creation implemented in create_junctions
    junctions_u$Top <- add_junction_heights(junctions_u, dtm)

    message("7. final junctions and conduits are created")

    #### 8. Add outfall connection of initial network ####
    conduits_sf$InOffset <- 0
    conduits_sf$OutOffset <- 0
    link_outfall <- conduits_sf[which(conduits_sf$ToNode %in% outfalls$Name),c("ToNode", "FromNode", "InOffset", "OutOffset", "L1")]

    for(n in 1:length(link_outfall$FromNode)){
      node_index <- which(junctions$Name == link_outfall$FromNode[n])
      dist <- sf::st_distance(junctions[node_index,], junctions_u)
      dist <- as.numeric(dist)
      link_outfall$FromNode[n] <- junctions_u$Name[dist == 0]
    }

    link_outfall$L1 <- seq(max(conduits_u$L1)+1, max(conduits_u$L1) + nrow(link_outfall), by = 1)
    conduits_u <- rbind(conduits_u, link_outfall)

    message("8. initial connections to outfalls are added to final network")

    #### Add name to conduits: ####
    conduits_u$Name <- paste0("C_", 1:nrow(conduits_u))

    #### Tag junctions: ####
    junctions_u <- error_checks(junctions = junctions_u, conduits_sf = conduits_u, outfalls = outfalls)

    # plot tags...
    plot_tagged_junctions(junctions_u, outfalls, conduits_u, main = "add tags to junctions")

    #### 9. Add conduits to drain remaining sinks ####
    if(any(junctions_u$tag == "sink")){
      count <- 1
      for(sink in junctions_u$Name[junctions_u$tag == "sink"]){
        # add a new conduit to ensure the permeability of the network:
        if(short_cut_sinks){
          to_sink <- conduits_u$FromNode[conduits_u$ToNode == sink]

          distances_sink <- sf::st_distance(junctions_u[junctions_u$Name == sink,], junctions_u)
          distances_sink <- as.numeric(distances_sink)
          distances_sink <- data.frame(distances = distances_sink, Name = junctions_u$Name, Top = junctions_u$Top, stringsAsFactors = F)

          distances_sink <- distances_sink[distances_sink$distance > 0 & !distances_sink$Name %in% to_sink &
                                             distances_sink$Top < junctions_u$Top[junctions_u$Name == sink],]

          name_nearest_lower_junction <-
            distances_sink$Name[distances_sink$distances == min(distances_sink$distances)]

          new_conduit <- rbind(junctions_u[junctions_u$Name == name_nearest_lower_junction,], junctions_u[junctions_u$Name == sink,])
          new_conduit <- dplyr::summarize(new_conduit)
          new_conduit <- sf::st_cast(new_conduit,"LINESTRING")

          new_conduit$ToNode <- name_nearest_lower_junction
        }

        # drain to the nearest outfall:
        if(direct_drainage_sinks){

          name_nearest_outfall <- sf::st_distance(junctions_u[junctions_u$Name == sink,], outfalls)
          name_nearest_outfall <-  as.numeric(name_nearest_outfall)
          name_nearest_outfall <- which.min(name_nearest_outfall)
          name_nearest_outfall <- outfalls$Name[name_nearest_outfall]

          # swmm does not allow two links at one outfall; therefore, add another outfall:
          new_outfall_name <- paste0("n_",name_nearest_outfall,"_",count)
          outfalls <- rbind(outfalls, outfalls[outfalls$Name == name_nearest_outfall,])
          outfalls[nrow(outfalls),"Name"] <- new_outfall_name
          outfalls[nrow(outfalls),"geometry"] <- outfalls[nrow(outfalls),"geometry"] + stats::runif(1, min = 0, max = 2)

          new_conduit <- rbind(outfalls[nrow(outfalls),"Name"], junctions_u[junctions_u$Name == sink,"Name"])
          new_conduit <- dplyr::summarize(new_conduit, .groups = "drop_last")
          new_conduit <- sf::st_cast(new_conduit, "LINESTRING")

          new_conduit$ToNode <- new_outfall_name

          count <- count + 1
        }

        new_conduit$Name <- paste0("nC_",nrow(conduits_u))
        new_conduit$FromNode <- sink
        new_conduit$L1 <- nrow(conduits_u) + 1
        new_conduit$InOffset <- 0
        new_conduit$OutOffset <- 0

        conduits_u <- rbind(new_conduit, conduits_u)

      }
    }

    message("9. Add short cuts to remaining local sinks")

    # plot tags...
    junctions_u <- error_checks(junctions = junctions_u, conduits_sf = conduits_u, outfalls = outfalls)
    plot_tagged_junctions(junctions_u, outfalls, conduits_u, main = "final error tags")
    plot_short_cuts(junctions_u, conduits_u, outfalls, "red")

    #### 10. Breaks in network (optional): ####
    message("10. optional breaks in network to increase linearity:")
    if(any(break_closed_loops, break_loops, breaks_at_hills, delete_disconnected)){
      list_breaks <- breaks_in_network(break_closed_loops, break_loops, breaks_at_hills, delete_disconnected, junctions = junctions_u, outfalls = outfalls, conduits_sf = conduits_u)
      junctions_u <- list_breaks$junctions
      conduits_u <- list_breaks$conduits
      outfalls <- list_breaks$outfalls
    }


    #### 11. Check consistency of network from outfall to start or hill ####
    # plot paths:
    for (i in 1:nrow(outfalls)){
      plot_path_to_outfall(junctions_u, outfalls, oid = i, conduits_u, main = paste0("path to outfall ", i), col_arrow = "red")
    }

    message("11. outfall catchments are plotted")

    #### 12. Slope and junction depth adjustment within the given range ####
    conduits_u$Length <- sf::st_length(conduits_u)
    conduits_u$Length <- as.numeric(conduits_u$Length)
    conduits_u$InOffset <- 0
    conduits_u$OutOffset <- 0

    # calculate initial slopes and plot them...
    conduits_u$Slope <- calculateSlope(junctions = junctions_u, conduits_sf = conduits_u, outfall = outfalls)
    graphics::boxplot(conduits_u$Slope, main = "conduit slopes before adjustment")
    graphics::boxplot(junctions_u$Top - junctions_u$Bottom, main = "junction depths before slope adjustment")
    plot_slopes(junctions = junctions_u, conduits_sf = conduits_u, min_slope = min_slope, max_slope = max_slope, outfalls = outfalls, main = "slopes before adjustment")

    # correct slopes and junction depths in the range min/max slope and junction depth
    list_slope_adjustments <- adjustSlopesAndJuncDepths(conduits = conduits_u, junctions = junctions_u, outfall = outfalls,
                                           min_slope = min_slope, max_slope = max_slope,
                                           min_junc_depth = min_junc_depth, max_junc_depth = max_junc_depth,
                                           mean_junc_depth = mean_junc_depth)
    outfalls <- list_slope_adjustments$outfall
    junctions_u <- list_slope_adjustments$junctions
    conduits_u <- list_slope_adjustments$conduits

    # calculate final slopes ..
    conduits_u$Slope <- calculateSlope(junctions = junctions_u, conduits_sf = conduits_u, outfall = outfalls)
    # plot slopes ...
    graphics::boxplot(conduits_u$Slope, main = "conduit slopes after adjustment")
    graphics::boxplot(junctions_u$Top - junctions_u$Bottom, main = "junction depths after slope adjustment")
    plot_slopes(junctions = junctions_u, conduits_sf = conduits_u, min_slope = min_slope, max_slope = max_slope, outfalls = outfalls, main = "slopes after adjustment")

    message("12. Slopes and junction depths are adjusted within the given range (min_junc_depth, max_junc_depth, min_slope, max_slope) ")

    #### 13. add remaining default parameters: ####
    conduits_u$Shape <- "CIRCULAR"
    conduits_u$Geom1 <- 0.3
    conduits_u$Roughness <- 0.018

    message("13. missing default parameters are added")

    #### return ####
    return(list(junctions = junctions_u, conduits = conduits_u, outfalls = outfalls))

  }


}
