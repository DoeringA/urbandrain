#' Create a SWMM model including a drainage network and recharging surfaces.
#'
#' Fully automated generation of a SWMM model including the drainage network and recharging surfaces based on open source spatial data. Requires street polylines, a digital terrain model, outfall point data, land cover data and some user defined parameters as input. First creates a drainage network following the location of the streets and the topography, then adds recharging surfaces and finally sized the pipe diameters while running a complete SWMM model. Stores SWMM input file and shp-files containing the drainage netwrok information in a defined directory. Returns a list of the sf objects junctions, conduits and outalls that can be converted to a SWMM *.inp file using swmmr functions afterwards.
#'
#' @param streets Street polyline data of class sf. On basis of the street polylines the drainage network is setup. Street polyline data can be extracted for example from openstreetmaps.
#' @param dtm Digital terrain information either of class RasterLayer or as a triangulated grid of class sf. Surface height data is required to add heights to junctions and to define the flow directions in the conduits subsequently.
#' @param outfalls Outfall point data of class sf. The network drains towards these points.
#' @param crs_default A proj4string or EPSG Code to specify the default coordinate reference system.
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
#' @param boundary_polygon A polygon file of class sf defining the model borders.
#' @param landuse_sf A polygon file of class sf including landuse informations.
#' @param landuse_classes A data.frame including percentages of imperviousness of the different landuse classes. See vignette.
#' @param path_timeseries Path to a raindata file in dat format including rain data to size the pipe diameter.
#' @param path_options Path to an option file specifying swmm options. See \link[swmmr]{shp_to_inp}
#' @param infiltration A data.frame specifying infiltration information for soil classes. Currently, homogenous soil for the model area is assumed.
#' @param path_out Directory to store the swmm input file and the model data in shp format to.
#' @return A list of sf objects junctions, conduits, outfalls and subcatchments.
#' @export
#' @rdname create_drainage_network

create_swmm_model <- function(streets, dtm, outfalls, crs_default, buffer, snap_dist, epsilon, lim,  min_junc_depth, mean_junc_depth, max_junc_depth,  min_slope, max_slope, ds, stepwise = T, break_closed_loops = F, delete_disconnected = F, breaks_at_hills = F, break_loops = F,  short_cut_sinks = F, direct_drainage_sinks = T, boundary_polygon, landuse_sf, landuse_classes, path_timeseries, path_options, infiltration, path_out){


  #... create drainage network
  network_list <- create_drainage_network(streets, dtm, outfalls, crs_default, buffer, snap_dist, epsilon, lim,  min_junc_depth, mean_junc_depth, max_junc_depth,  min_slope, max_slope, ds, stepwise, break_closed_loops,  delete_disconnected, breaks_at_hills, break_loops,  short_cut_sinks, direct_drainage_sinks)

  junctions <- network_list$junctions
  conduits <- network_list$conduits
  outfalls <- network_list$outfalls

  # ... create SWMM subcatchments
  if(exists("landuse_sf")){
    landuse_sf <- sf::st_transform(landuse_sf, crs_default)
    boundary_polygon <- sf::st_transform(boundary_polygon, crs_default)
    # create subcatchments with thiessen polygons and landuse information:
    subcatchments <- thiessenpolygons(junctions, boundary_polygon)
    subcatchments <- thiessen_to_subcatchments(landuse_sf = landuse_sf, landuse_classes = landuse_classes, thiessen_polygons = subcatchments)
  }

  # # ... delete all files except options and rain:
  # outfiles <- list.files(path_out, full.names = TRUE)[!list.files(path_out, full.names = TRUE) %in% list(file.path(path_out, "default_rain.dat"), file.path(path_out, "options.txt"))]
  # do.call(file.remove, list(outfiles))

  # # ... save sf to shp
  # sf::st_write(junctions, file.path(path_out, "nodes_artificial_SWMM_format.shp"))
  # sf::st_write(conduits, file.path(path_out, "links_artificial_SWMM_format.shp"))
  # sf::st_write(outfalls, file.path(path_out, "outfall_artificial_SWMM_format.shp"))
  # sf::st_write(subcatchments, file.path(path_out, "subcatchments_artificial_SWMM_format.shp"))

  # ... convert sf to inp
  inp <- swmmr::sf_to_inp(
    outfall_sf = outfalls,
    line_sf = conduits,
    point_sf = junctions,
    polygon_sf = subcatchments,
    path_timeseries = path_timeseries,
    path_options = path_options,
    infiltration = infiltration
  )

  # ... save inp
  swmmr::write_inp(inp, file.path(path_out, "artificial_SWMM_model.inp"))

  # ... run swmm
  swmmr::run_swmm(file.path(path_out, "artificial_SWMM_model.inp"))

  # ... pipe sizing algorithm
  conduits_sized <- pipeSizingAlgorithm(inp, path_out, conduits, junctions, outfalls, target = "conduit_surcharge")

  # # ... save shp of sized pipes
  # sf::st_write(conduits_sized, file.path(path_out, "final_links_artificial_SWMM_format.shp"))

  # return SWMM network data
  return(list(junctions = junctions, conduits = conduits_sized, outfalls = outfalls, subcatchments = subcatchments, inp = inp))
}

