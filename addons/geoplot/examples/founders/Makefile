gpfiles: founders.inp founders.geojson founders.csv
	gretlcli -t founders.inp

gnuplot: gpfiles
	gnuplot -persist GDPpc.gp
clean:
	rm -f *~ session.inp string_table.txt
	rm -f GDPpc.* tmp.dat
