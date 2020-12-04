#' plot of drainage network during creation process
#' @keywords internal
plot_drainage_network <- function(coords, outfalls, conduits_sf, main, col_arrow){
  # initiate plot
  graphics::plot(coords$X, coords$Y, xlab = "X", ylab = "Y", main = main)
  coords_O <- as.data.frame(sf::st_coordinates(outfalls))
  coords_O$Name <- outfalls$Name
  graphics::points(coords_O$X, coords_O$Y, col = "green", pch = 19)


  for(i in 1:nrow(conduits_sf)){
    p0 <- conduits_sf$FromNode[i]
    p1 <- conduits_sf$ToNode[i]
    # plot arrow:
    if(all(list(p0,p1) %in% coords$Name)){
      graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
             x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
             col = col_arrow, length = 0.1)
    }
    if(all(list(p0,p1) %in% coords_O$Name)){
      graphics::arrows(x0 = coords_O$X[coords_O$Name == p0], y0 = coords_O$Y[coords_O$Name == p0],
             x1 = coords_O$X[coords_O$Name == p1], y1 = coords_O$Y[coords_O$Name == p1],
             col = col_arrow, length = 0.1)
    }
    if(p0 %in% coords$Name & p1 %in% coords_O$Name){
      graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
             x1 = coords_O$X[coords_O$Name == p1], y1 = coords_O$Y[coords_O$Name == p1],
             col = col_arrow, length = 0.1)
    }
    if(p0 %in% coords_O$Name & p1 %in% coords$Name){
      graphics::arrows(x0 = coords_O$X[coords_O$Name == p0], y0 = coords_O$Y[coords_O$Name == p0],
             x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
             col = col_arrow, length = 0.1)
    }

  }
}

#' plot tagged junctions
#' @keywords internal
plot_tagged_junctions <- function(junctions, outfalls, conduits_sf, main){
  coords <- as.data.frame(sf::st_coordinates(junctions))
  coords$Name <- junctions$Name
  coords$tag <- junctions$tag

  dic <- data.frame(tag = c("hill", "sink", "outfall_artificial", "loop", "normal", "start",  "crossing"),
                    pch = c(17,15,19,19,1,19,19), col = c("red", "red", "red", "lightblue", "black", "black", "blue"),
                    cex = c(1.3,1.3,1.3,1,1,1,1),
                    stringsAsFactors = F)

  graphics::par(mar=c(5.1, 4.1, 4.1, 9), xpd=TRUE)
  plot_drainage_network(coords, outfalls, conduits_sf, main = main, col_arrow = "black")
  for(tag in unique(coords$tag)){
    graphics::points(coords$X[coords$tag == tag], coords$Y[coords$tag == tag], pch = dic$pch[dic$tag == tag],
           col = dic$col[dic$tag == tag], cex = dic$cex[dic$tag == tag])
  }
  graphics::legend("topright",inset = c(-0.4,0), pch = c(19, dic$pch), col = c("green", dic$col), c("outfalls", dic$tag))
}

#' plot catchment of outfalls
#' @keywords internal
plot_path_to_outfall <- function(junctions_u, outfalls, oid, conduits_u, main, col_arrow){

  # preparation default plot:
  coords <- as.data.frame(sf::st_coordinates(junctions_u))
  coords$Name <- junctions_u$Name

  # default plot:
  plot_drainage_network(coords, outfalls, conduits_u, main = main, col_arrow = "black")

  # colour selected outfall:
  coords_O <- as.data.frame(sf::st_coordinates(outfalls[oid,]))
  coords_O$Name <- outfalls$Name[oid]
  graphics::points(coords_O$X, coords_O$Y, col = "red", pch = 19)

  # paths to outfall
  list_O <- track_outfall_to_start_connectivity(outfalls = outfalls[oid,], junctions = junctions_u, conduits_sf = conduits_u)

  # add coloured arrows to show catchment of a defined outfall
  for(t in 1:length(list_O)){
    for(j in 1:length(list_O[[t]])){
      to_O <- c(outfalls$Name[oid], names(list_O[[t]][[j]]))

      for(i in 1:(length(to_O)-1)){
        p0 <- to_O[i+1]
        p1 <- to_O[i]
        # plot arrows:
        if(all(list(p0,p1) %in% coords$Name)){
          graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
                 x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
                 col = col_arrow, length = 0.1)
        }
        if(all(list(p0,p1) %in% coords_O$Name)){
          graphics::arrows(x0 = coords_O$X[coords_O$Name == p0], y0 = coords_O$Y[coords_O$Name == p0],
                 x1 = coords_O$X[coords_O$Name == p1], y1 = coords_O$Y[coords_O$Name == p1],
                 col = col_arrow, length = 0.1)
        }
        if(p0 %in% coords$Name & p1 %in% coords_O$Name){
          graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
                 x1 = coords_O$X[coords_O$Name == p1], y1 = coords_O$Y[coords_O$Name == p1],
                 col = col_arrow, length = 0.1)
        }
        if(p0 %in% coords_O$Name & p1 %in% coords$Name){
          graphics::arrows(x0 = coords_O$X[coords_O$Name == p0], y0 = coords_O$Y[coords_O$Name == p0],
                 x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
                 col = col_arrow, length = 0.1)
        }

      }

    }
  }

}

#' cross section plot when adjusting pipe slopes and junction depths
#' @keywords internal
plot_cross_section_junctions <- function(junctions, conduits, connections){
  cross_section_data <- data.frame(Name = connections, Length = 0, Offset = 0, stringsAsFactors = F)

  for(i in 1:length(connections)){
    name <- connections[i]
    cross_section_data$Top[i] <- junctions$Top[junctions$Name == name]
    cross_section_data$Bottom[i] <- junctions$Bottom[junctions$Name == name]
    cross_section_data$max_bottom[i] <- junctions$max_bottom[junctions$Name == name]
    cross_section_data$min_bottom[i] <- junctions$min_bottom[junctions$Name == name]
  }

  for(i in 1:(length(connections)-1)){
    length_part <- conduits$Length[conduits$FromNode == connections[i+1] & conduits$ToNode == connections[i]]
    cross_section_data$Length[i+1] <- cross_section_data$Length[i] + length_part
    cross_section_data$Offset[i] <- conduits$OutOffset[conduits$FromNode == connections[i+1] & conduits$ToNode == connections[i]]
  }


  cross_section_data$Bottom_pipe <- cross_section_data$Bottom + cross_section_data$Offset


  graphics::plot(cross_section_data$Length, cross_section_data$Bottom_pipe,
       ylim = c(min(cross_section_data$Bottom)-5, max(cross_section_data$Top)+1),
       xlab = "Length [m]",
       ylab = "Height [m a.s.l.]",
       type = "b", pch = 19)

  graphics::text(cross_section_data$Length, cross_section_data$Top + 0.3, cross_section_data$Name, cex = 0.5, col = "grey")

  graphics::rect(xleft = cross_section_data$Length - 1, ybottom = cross_section_data$Bottom,
       xright = cross_section_data$Length + 1, ytop = cross_section_data$Top,
       col = "lightgrey", border = NA)

  graphics::points(cross_section_data$Length, cross_section_data$Top, type = "b", pch = 19)
  graphics::points(cross_section_data$Length, cross_section_data$Bottom, pch = 19, col = "grey")
  graphics::points(cross_section_data$Length, cross_section_data$Bottom_pipe, type = "b", pch = 19)

  graphics::points(cross_section_data$Length, cross_section_data$max_bottom, col = "red", pch = 4)
  graphics::points(cross_section_data$Length, cross_section_data$min_bottom, col = "red", pch = 4)

  graphics::legend("bottomright", c("top and bottom ", "possible range of bottom", "offset"), pch = c(19, 4, 19), col = c("black", "red", "grey"))

}

#' plot slopes
#' @keywords internal
plot_slopes <- function(junctions, outfalls, conduits_sf, min_slope, max_slope, main, roundValue = 3){

  conduits_sf$Slope <- round(conduits_sf$Slope, roundValue)

  conduits_sf$col_arrow <- "black"
  conduits_sf$col_arrow[conduits_sf$Slope < min_slope] <- "red"
  conduits_sf$col_arrow[conduits_sf$Slope > max_slope] <- "purple"

  # coordinates from juncions:
  coords <- as.data.frame(sf::st_coordinates(junctions))
  coords$Name <- junctions$Name

  # initiate plot:
  graphics::par(mar=c(5.1, 4.1, 4.1, 9), xpd=TRUE)

  graphics::plot(coords$X, coords$Y, xlab = "X", ylab = "Y", main = main)
  coords_O <- as.data.frame(sf::st_coordinates(outfalls))
  coords_O$Name <- outfalls$Name
  graphics::points(coords_O$X, coords_O$Y, col = "black", pch = 19)

  for(i in 1:nrow(conduits_sf)){
    p0 <- conduits_sf$FromNode[i]
    p1 <- conduits_sf$ToNode[i]
    # plot arrow:
    if(all(list(p0,p1) %in% coords$Name)){
      graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
             x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
             col = conduits_sf$col_arrow[i], length = 0.1)
    }
    if(all(list(p0,p1) %in% coords_O$Name)){
      graphics::arrows(x0 = coords_O$X[coords_O$Name == p0], y0 = coords_O$Y[coords_O$Name == p0],
             x1 = coords_O$X[coords_O$Name == p1], y1 = coords_O$Y[coords_O$Name == p1],
             col = conduits_sf$col_arrow[i], length = 0.1)
    }
    if(p0 %in% coords$Name & p1 %in% coords_O$Name){
      graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
             x1 = coords_O$X[coords_O$Name == p1], y1 = coords_O$Y[coords_O$Name == p1],
             col = conduits_sf$col_arrow[i], length = 0.1)
    }
    if(p0 %in% coords_O$Name & p1 %in% coords$Name){
      graphics::arrows(x0 = coords_O$X[coords_O$Name == p0], y0 = coords_O$Y[coords_O$Name == p0],
             x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
             col = conduits_sf$col_arrow[i], length = 0.1)
    }
  }

  # add plot legend:
  graphics::legend("topright", inset = c(-0.4,0), legend = c("slope < min_slope", "slope > max_slope"), col = c("red", "purple"), lwd = 1)

}


#' plot conduits that are not connected to an outfall
#' @keywords internal
plot_disconnected_conduits <- function(junctions, conduits, outfalls, col_arrow, main = "conduits disconnected from outfalls"){

  # preparation default plot:
  coords <- as.data.frame(sf::st_coordinates(junctions))
  coords$Name <- junctions$Name

  # extract outfall coordinates:
  coords_O <- as.data.frame(sf::st_coordinates(outfalls))
  coords_O$Name <- outfalls$Name


  # default plot:
  plot_drainage_network(coords, outfalls, conduits, main = main, col_arrow = "black")

  # paths to outfall
  network_paths <- track_outfall_to_start_connectivity(outfalls, junctions, conduits)

  all_conduits_connected <- NULL
  for(j in 1:length(network_paths)){
    for(i in 1:length(network_paths[[j]])){
      names_junc <- names(network_paths[[j]][[i]])

      conduits_connected <- NULL
      for(n in 2:length(names_junc)){
        to <- names_junc[n-1]
        from <- names_junc[n]
        c_name <- conduits$Name[conduits$FromNode == from & conduits$ToNode == to]
        if(length(c_name) > 1){
          message(paste("junctions are connected more than once:", paste(c_name, collapse = " ")))
          c_name <- c_name[1]
        }
        conduits_connected[n-1] <- c_name

      }


      if(is.null(all_conduits_connected)){
        all_conduits_connected <- conduits_connected
      }else{
        all_conduits_connected <- c(all_conduits_connected, conduits_connected)
      }
    }
  }

  conduits_disconnected <- conduits$Name[!(conduits$Name %in% unique(all_conduits_connected) | conduits$ToNode %in% outfalls$Name)]

  if(length(conduits_disconnected)>0){
    # add coloured arrows to mark disconnected conduits
    for(i in 1:length(conduits_disconnected)){
      p0 <- conduits$FromNode[conduits$Name == conduits_disconnected[i]]
      p1 <- conduits$ToNode[conduits$Name == conduits_disconnected[i]]
      # plot arrows:
      if(all(list(p0,p1) %in% coords$Name)){
        graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
               x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
               col = col_arrow, length = 0.1)
      }
      if(all(list(p0,p1) %in% coords_O$Name)){
        graphics::arrows(x0 = coords_O$X[coords_O$Name == p0], y0 = coords_O$Y[coords_O$Name == p0],
               x1 = coords_O$X[coords_O$Name == p1], y1 = coords_O$Y[coords_O$Name == p1],
               col = col_arrow, length = 0.1)
      }
      if(p0 %in% coords$Name & p1 %in% coords_O$Name){
        graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
               x1 = coords_O$X[coords_O$Name == p1], y1 = coords_O$Y[coords_O$Name == p1],
               col = col_arrow, length = 0.1)
      }
      if(p0 %in% coords_O$Name & p1 %in% coords$Name){
        graphics::arrows(x0 = coords_O$X[coords_O$Name == p0], y0 = coords_O$Y[coords_O$Name == p0],
               x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
               col = col_arrow, length = 0.1)
      }

    }
  }



}

#' plot short cuts
#' @keywords internal
plot_short_cuts <- function(junctions, conduits, outfalls, col_arrow){
  # preparation default plot:
  coords <- as.data.frame(sf::st_coordinates(junctions))
  coords$Name <- junctions$Name

  # extract outfall coordinates:
  coords_O <- as.data.frame(sf::st_coordinates(outfalls))
  coords_O$Name <- outfalls$Name


  # default plot:
  plot_drainage_network(coords, outfalls, conduits, main = "conduits added to drain sinks", col_arrow = "black")

  # added short cuts:
  added_conduits <- conduits$Name[grep("n", conduits$Name)]

  if(length(added_conduits > 0)){
    # add coloured arrows to mark disconnected conduits
    for(i in 1:length(added_conduits)){
      p0 <- conduits$FromNode[conduits$Name == added_conduits[i]]
      p1 <- conduits$ToNode[conduits$Name == added_conduits[i]]
      # plot arrows:
      if(all(list(p0,p1) %in% coords$Name)){
        graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
               x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
               col = col_arrow, length = 0.1)
      }
      if(all(list(p0,p1) %in% coords_O$Name)){
        graphics::arrows(x0 = coords_O$X[coords_O$Name == p0], y0 = coords_O$Y[coords_O$Name == p0],
               x1 = coords_O$X[coords_O$Name == p1], y1 = coords_O$Y[coords_O$Name == p1],
               col = col_arrow, length = 0.1)
      }
      if(p0 %in% coords$Name & p1 %in% coords_O$Name){
        graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
               x1 = coords_O$X[coords_O$Name == p1], y1 = coords_O$Y[coords_O$Name == p1],
               col = col_arrow, length = 0.1)
      }
      if(p0 %in% coords_O$Name & p1 %in% coords$Name){
        graphics::arrows(x0 = coords_O$X[coords_O$Name == p0], y0 = coords_O$Y[coords_O$Name == p0],
               x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
               col = col_arrow, length = 0.1)
      }

    }

  }


}

#' plot pipe diameters
#' @keywords internal
plot_drainage_network_pipe_diameters <- function(junctions, outfalls, conduits_sf, main, DN_m){

  # scale arrow line thickness according to diameter:
  dic_scaling <- data.frame(value = DN_m[order(DN_m, decreasing = T)], stringsAsFactors = F)
  dic_scaling$col_arrow <- viridis::viridis(nrow(dic_scaling))
  dic_scaling$lwd <- round(seq(10, 1, length.out = nrow(dic_scaling)),1) # maximum line width 10 ?

  # initialize columns:
  conduits_sf$lwd <- NA
  conduits_sf$col_arrow <- NA

  # add colors to conduits:
  for(i in 1:nrow(conduits_sf)){
    conduits_sf$lwd[i] <- dic_scaling$lwd[as.character(dic_scaling$value) == as.character(conduits_sf$Geom1[i])]
    conduits_sf$col_arrow[i] <- dic_scaling$col_arrow[as.character(dic_scaling$value) == as.character(conduits_sf$Geom1[i])]
  }

  # coordinates from juncions:
  coords <- as.data.frame(sf::st_coordinates(junctions))
  coords$Name <- junctions$Name

  # initiate plot:
  graphics::par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

  graphics::plot(coords$X, coords$Y, xlab = "X", ylab = "Y", main = main)
  coords_O <- as.data.frame(sf::st_coordinates(outfalls))
  coords_O$Name <- outfalls$Name
  graphics::points(coords_O$X, coords_O$Y, col = "black", pch = 19)

  for(i in 1:nrow(conduits_sf)){
    p0 <- conduits_sf$FromNode[i]
    p1 <- conduits_sf$ToNode[i]
    # plot arrow:
    if(all(list(p0,p1) %in% coords$Name)){
      graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
             x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
             col = conduits_sf$col_arrow[i], length = 0.1, lwd = conduits_sf$lwd[i])
    }
    if(all(list(p0,p1) %in% coords_O$Name)){
      graphics::arrows(x0 = coords_O$X[coords_O$Name == p0], y0 = coords_O$Y[coords_O$Name == p0],
             x1 = coords_O$X[coords_O$Name == p1], y1 = coords_O$Y[coords_O$Name == p1],
             col = conduits_sf$col_arrow[i], length = 0.1, lwd = conduits_sf$lwd[i])
    }
    if(p0 %in% coords$Name & p1 %in% coords_O$Name){
      graphics::arrows(x0 = coords$X[coords$Name == p0], y0 = coords$Y[coords$Name == p0],
             x1 = coords_O$X[coords_O$Name == p1], y1 = coords_O$Y[coords_O$Name == p1],
             col = conduits_sf$col_arrow[i], length = 0.1, lwd = conduits_sf$lwd[i])
    }
    if(p0 %in% coords_O$Name & p1 %in% coords$Name){
      graphics::arrows(x0 = coords_O$X[coords_O$Name == p0], y0 = coords_O$Y[coords_O$Name == p0],
             x1 = coords$X[coords$Name == p1], y1 = coords$Y[coords$Name == p1],
             col = conduits_sf$col_arrow[i], length = 0.1, lwd = conduits_sf$lwd[i])
    }
  }

  # add plot legend:
  content_legend <- dic_scaling[dic_scaling$value %in% unique(conduits_sf$Geom1),]

  graphics::legend("topright", inset = c(-0.3,0), legend = content_legend$value, lwd = content_legend$lwd, col = content_legend$col_arrow)

}
