CFLAGS = -g -Wall -I/opt/gnome/include/gnome-xml `gtk-config --cflags`
LDFLAGS = -L/opt/gnome/lib -lxml `gtk-config --libs`

all: gnumeric_import


