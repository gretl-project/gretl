GnomeUIInfo gnome_helpmenu[] = {
    {GNOME_APP_UI_ITEM, 
     N_("About"), N_("About this program"),
     about_cb, NULL, NULL, 
     GNOME_APP_PIXMAP_STOCK, GNOME_STOCK_MENU_ABOUT,
     0, 0, NULL},
    GNOMEUIINFO_SEPARATOR,
    GNOMEUIINFO_HELP("gretl"),
    GNOMEUIINFO_END
};

