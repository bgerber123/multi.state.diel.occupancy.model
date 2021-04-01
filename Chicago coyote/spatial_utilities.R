check_path <- function(path){
  if( file.exists( path ) ){
    cat(
      paste(
        cli::col_green(
          cli::symbol$tick
        ),
        "Raster layer exists")
    )
  }else{
    # The error (warning to provide to the user.)
    err <- paste0(
      "File does not exist.\n",
      "  You have either entered the wrong path or\n",
      "  you need to download the appropriate raster."
    )
    cat(
      paste(
        cli::col_red(
          cli::symbol$cross
        ),
        err
      )
    )
  }
}

# a hidden function for some command line reporting
.cli_post <- function(my_text, pass = FALSE){
  if(pass){
    cat(
      paste(
        cli::col_green(
          cli::symbol$tick
        ),
        "\n"
      )
    )
    
  } else {
    cat(
      paste(
        cli::symbol$bullet,
        my_text,
        " "
      )
    )
  }
}

# extract land use landcover data

extract_raster_prop <- function(
  my_points,
  my_buffer = NULL,
  my_raster_data,
  lulc_cats = NULL,
  point_names = NULL
){
  # This function has 5 arguments:
  #  my_points (sf object): The coordinates  of the sites you want to collect
  #                           spatial data from. The column name for the
  #                           names of the sites must be 'LocationName'.
  #  my_buffer (numeric):   The radius of a buffer around each point form which
  #                           to extract cell values. If the distance between 
  #                           the sampling point and the center of a cell is 
  #                           less than or equal to the buffer, the cell is 
  #                           included. The buffer can be specified as a single
  #                           value, or a vector of the length of the number of
  #                           points. If the data are not projected (e.g.,
  #                           latitude/longitude), unit should be in meters.
  #                           Otherwise it should be in map units, which is also
  #                           typically meters.
  #  my_raster_data (raster): A raster object that you want to extract data 
  #                             from.
  #  lulc_cats (numeric or list, optional): Land-use land-cover categories. 
  #                                 If this is a numeric then the data will
  #                                 be queried to return just those specific
  #                                 lulc categories (e.g., lulc_cats = c(1,2))
  #                                 will return the first two categories. If
  #                                 this is a list then each element must be a 
  #                                 numeric vector. If the numeric vector is
  #                                 > 1 for a specific list element then those
  #                                 categories will be summed together. This
  #                                 may be done to combine specific similar
  #                                 category types (e.g., combine 'building',
  #                                 'road', and 'other paved surfaces') to 
  #                                 generate an impervious surface category.
  #                                 The numeric or the list can be named. If
  #                                 they are then those names will be added
  #                                 to the returned object. If NULL, then all
  #                                 categories will be returned.
  #  point_names (character or factor): The names of the point locations, 
  #                                 this is generally a column in the 
  #                                 my_points object (e.g., my_points$NAME)
  #                                 supplied to the first function argument.
  # start command line reporting
  cat(
    paste0(
      cli::rule(line = 2),
      "\n\n"
    )
  )
  # Step 1. reproject the point data to match raster projection
  .cli_post("Reprojecting my_points to map projection:")
  
  points_RP <- sf::st_transform(
    my_points, 
    sf::st_crs(my_raster_data)
  )
  
  .cli_post(pass = TRUE)
  
  # Step 2.
  # use the raster::extract function to extract the mean value of the raster data 
  # within a particular buffer around each site.
  # We need a sub-function to calculate the proportion of each category
  spatial_summary <- function(
    x,
    ncats = my_raster_data@data@max,
    ...
  ){
    return(
      prop.table(
        tabulate(x, ncats)
      )
    )
  }
  
  .cli_post("Extracting spatial data:")
  prop_extract <- raster::extract(
    my_raster_data,
    points_RP,
    fun=spatial_summary,
    buffer= my_buffer,
    na.rm = TRUE
  )
  .cli_post(pass = TRUE)
  
  if(is.numeric(prop_extract)){
  prop_extract <- matrix(
    prop_extract,
    ncol = my_raster_data@data@max,
    nrow = nrow(points_RP),
    byrow = TRUE
  )
  }
  
  #if(!is.null(point_names)){
  #row.names(prop_extract) <- point_names
  #}
  
  # if lulc_cats is a list
  if(is.list(lulc_cats)){    
    prop_extract <- apply(
      prop_extract,
      1,
      function(x){
        sapply(
          lulc_cats,
          function(y)
            ifelse(
              length(y) == 1,
              x[y],
              sum(x[y])
            )
        )
      }
    )
    prop_extract <- t(prop_extract)
  }
  
  # if it is a numeric
  if(is.numeric(lulc_cats)){
    prop_extract <- prop_extract[,lulc_cats]
    if(!is.null(names(lulc_cats))){
      colnames(prop_extract) <- names(lulc_cats)
    }
  }
  
  # create dataframe matching the sites with the extracted data
  .cli_post("Summarizing data:")
  df <- data.frame(
    LocationName = point_names,
    #  if(is.matrix(prop_extract)){
    #  row.names(prop_extract)
    #} else {
    #  names(prop_extract)
    #},
    prop_extract,
    row.names = NULL
  )
  .cli_post(pass = TRUE)
  # end command line reporting
  cat(
    paste0(
      "\n",
      cli::rule(line = 2),
      "\n\n"
    )
  )
  return(df[order(df$LocationName),])
}



# extract from a polygon

extract_polygon <- function(
  my_points,
  my_buffer,
  my_shp,
  layers = NULL
){
  ##############################################################################
  # This function has 4 arguments:
  #  my_points (sf object): The coordinates  of the sites you want to collect
  #                           spatial data from. The column name for the
  #                           names of the sites must be 'LocationName'.
  #  my_buffer (numeric):   The radius of a buffer around each point form which
  #                           to extract cell values. If the distance between 
  #                           the sampling point and the center of a cell is 
  #                           less than or equal to the buffer, the cell is 
  #                           included. The buffer can be specified as a single
  #                           value, or a vector of the length of the number of
  #                           points. If the data are not projected (e.g.,
  #                           latitude/longitude), unit should be in meters.
  #                           Otherwise it should be in map units, which is also
  #                           typically meters.
  #  my_shp (sf):           A sf object of a shape file you want to extract
  #                           data from.
  #  layers (character):    A character vector of column names from my_shp
  #                           to summarise. Optional. If NULL, all numeric
  #                           columns are summarised.
  ##############################################################################
  
  # If an object is supplied to layers, make sure it is a character.
  warn <- NULL
  if(!is.null(layers)){
    if(!is.character(layers)){
      stop("A non-character object has been supplied to layers.")
    }
    # check if all the layers exist, warn if not.
    if(!all(layers %in% colnames(my_shp))){
      warn <- paste("\nWarning: 'my_shp' did not have all layers specified in",
                    "'layers' argument.\nMissing elements in 'layers' were",
                    "removed'.")
      # drop the elements in layer that lack columns in my_shp
      layers <- layers[which(layers %in% colnames(my_shp))]
      if(length(layers) == 0){
        layers <- NULL
      }
    }
  }
  # STEP 1
  # Reproject the point data to match projection of population layer
  .cli_post("Reprojecting my_points to map projection:")
  
  points_RP <- sf::st_transform(
    my_points,
    sf::st_crs(my_shp)
  )
  
  .cli_post(pass = TRUE)
  # STEP 2 Create a buffer around each site
  .cli_post("Extracting spatial data:")
  
  points_buffer <- sf::st_buffer(
    points_RP,
    dist = my_buffer
  )
  
  # Extract population layer features that are contained within each buffers
  shp_intersection <- sf::st_intersection(
    points_buffer,
    my_shp
  )
  
  # Summarise the data if layers is not null
  if(!is.null(layers)){
    summary_data <- shp_intersection %>% 
      # remove spatial geometry so you are left with a data frame
      sf::st_set_geometry(NULL) %>% 
      # group by Location Name
      group_by(LocationName) %>% 
      # sum all the intersected pieces to get total housing units in each buffer
      dplyr::summarise_at(.vars = layers, .funs = sum) %>% 
      # divide by area of buffer, converting to km^2
      #  spatial data are the clumns we are summarizing
      #  sf::st_area(points_buffer) is the area we are buffering
      #  we modify the units 
      dplyr::mutate_if(is.numeric, 
                       function(spatial_data){
                         spatial_data/units::set_units(
                           sf::st_area(
                             points_buffer
                           ),
                           km^2
                         )
                       }
      ) %>% 
      arrange(LocationName)
    
  } else {
    summary_data <- shp_intersection %>% 
      # remove spatial geometry so you are left with a data frame
      sf::st_set_geometry(NULL) %>% 
      # group by Location Name
      group_by(LocationName) %>% 
      # sum all the intersected pieces to get total housing units in each buffer
      dplyr::summarise_if(is.numeric, .funs = sum) %>% 
      # divide by area of buffer, converting to km^2
      #  spatial data are the clumns we are summarizing
      #  sf::st_area(points_buffer) is the area we are buffering
      #  we modify the units 
      dplyr::mutate_if(is.numeric, 
                       function(spatial_data){
                         spatial_data/units::set_units(
                           sf::st_area(
                             points_buffer
                           ),
                           km^2
                         )
                       }
      )
    
  }
  
  .cli_post(pass = TRUE)
  if(!is.null(warn)){
    cat(
      cli::bg_black(
        cli::col_yellow(
          warn
        )
      )
    )
  }
  return(summary_data)
}


