# geoplot resources

# map files
founders.geojson     founder countries of the EU
switzerland.geojson  Switzerland, showing cantons
states.shp           US states and DC (shapefile)
us-states.geojson    US states and DC
us2020.geojson       US states and DC with added properties

# data files
founders.csv   population and GDP, EU founders (cross section, n=6)
statepop.gdt   US states and DC, population (cross section, n=51)

# scripts
edit_example.inp     shrinkage, placement, simplification of features[*]
founders.inp         plot GDP per capita, original EU countries
founders_mod.inp     add GDP data to geojson file[*]
seek_example.inp     find named features in map file
swiss-langs.inp      Swiss map: main languages per canton
us-2020.inp          Outcomes of 2020 US Presidential election[*]
us.inp               US map, random data
us-pop-density.inp   US map, showing population density[*]
us_states_json.inp   US map, polygons only
us_states_shp.inp    variant using shapefile input

# The scripts marked with [*] write out geojson files
