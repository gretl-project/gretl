<?php

$gretl_version = $_GET["gretl_version"];

$gretl_maj = strtok($gretl_version, ".");
$gretl_min = strtok(".");
$gretl_pl = strtok(".");

if ($dh = opendir(".")) {
    while (($file = readdir($dh)) !== false) {
        if (strstr($file, ".versions")) {
           $pkg = strstr($file, ".versions", TRUE);
           echo "found versions file for $pkg \n";
        }
    }
    closedir($dh);
}

?>

