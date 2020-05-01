all: founders usa

founders: founders.inp founders.geojson founders.csv
	gretlcli -t founders.inp

usa: us.inp
	gretlcli -t us.inp

blue: states.shp states.dbf states.shx us_states_shp.inp
	gretlcli -t us_states_shp.inp

outlines: 
	gretlcli -t us_states_json.inp

clean:
	rm -f string_table.txt *.gp *.dat *~
