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

$gretl_version = $_GET["gretl_version"];
$pkg = $_GET["pkg"];

$gretl_maj = strtok($gretl_version, ".");
$gretl_min = strtok(".");
$gretl_pl = strtok(".");

$version_data = $pkg . ".versions";
$fp = fopen($version_data, "r") or exit("Unable to open " . $version_data . "!");

$last_version = "";
$got_version = "";

while (!feof($fp)) {
    $pkg_version = trim(fgets($fp));
    if (strlen($pkg_version) > 0) {
	$vp1 = strtok($pkg_version, ".");
	$vp2 = strtok(".");
	$vp3 = strtok(".");
	if (too_new($gretl_maj, $gretl_min, $gretl_pl, $vp1, $vp2, $vp3)) {
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
    echo "pkgdir for ".$pkg." for gretl ".$gretl_version.": ".$got_version. "\n";
} else {
    echo "pkgdir for ".$pkg." for gretl ".$gretl_version.": none\n";
}

?>

