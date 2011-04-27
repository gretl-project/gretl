<?php
function too_new ($vg1, $vg2, $vg3, $vp1, $vp2, $vp3)
{
    if ($vp1 > $vg1) {
	return 1;
    } elseif ($vp1 == $vg1 && $vp2 > $vg2) {
	return 1;
    } elseif ($vp1 == $vg1 && $vp2 == $vg2 && $vp3 > $vg3) {
	return 1;
    } else {
	return 0;
    }
}

function get_latest ($file, $pkg, $vg1, $vg2, $vg3)
{
    $fp = @fopen($file, "r");
 
    if (@fp == false) {       
        return 0;
    }

    $last_version = "";
    $got_version = "";

    while (($line = fgets($fp, 32)) !== false) {
        $pkg_version = trim($line);
        if (strlen($pkg_version) > 0) {
	    $vp1 = strtok($pkg_version, ". ");
	    $vp2 = strtok(". ");
	    $vp3 = strtok(". ");
	    if (too_new($vg1, $vg2, $vg3, $vp1, $vp2, $vp3)) {
	        break;
	    } else {
	        $last_version = $pkg_version;
	    }
        }
    }

    fclose($fp);

    if (strlen($last_version) > 0) {
        $got_version = $last_version;
    }
    if (strlen($got_version)) {
        echo "$pkg $got_version\n";
        // $specfile = "$pkg.spec";
        // $spec = file_get_contents($specfile);
        // if ($spec !== false) {
        //     echo $spec;
        // }
        return 1;
    } 
    return 0;
}

$gretl_version = $_GET["gretl_version"];
$gretl_maj = strtok($gretl_version, ".");
$gretl_min = strtok(".");
$gretl_pl = strtok(".");

if ($dh = opendir(".")) {
    while (($file = readdir($dh)) !== false) {
        if (is_file($file) and strstr($file, ".versions")) {
           $pkgname = strstr($file, ".versions", TRUE);
           // echo "*** found versions file for $pkgname \n";
           get_latest($file, $pkgname, $gretl_maj, $gretl_min, $gretl_pl);
        }
    }
    closedir($dh);
}

?>

