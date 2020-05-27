all: founders usa blue outlines edit seek 

founders: founders.inp founders.geojson founders.csv
	gretlcli -t founders.inp

usa: us.inp
	gretlcli -t us.inp

blue: states.shp states.dbf states.shx us_states_shp.inp
	gretlcli -t us_states_shp.inp

outlines: 
	gretlcli -t us_states_json.inp

seek:
	gretlcli -t seek_example.inp 

edit:
	gretlcli -t edit_example.inp 

clean:
	rm -f string_table.txt
	rm -f *.plt *.dat *~

