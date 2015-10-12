<?php
function letter_to_int ($c)
{
    if ($c == 'a') return 0;
    if ($c == 'b') return 1;
    if ($c == 'c') return 2;
    if ($c == 'd') return 3;
    if ($c == 'e') return 4;
    if ($c == 'f') return 5;
    if ($c == 'g') return 6;
    if ($c == 'h') return 7;
    if ($c == 'i') return 8;
    if ($c == 'j') return 9;
    return 0;
}

$gretl_version = $_GET["gretl_version"];
$pkg = $_GET["pkg"];

$n = sscanf($gretl_version, "%d.%d.%d", $g1, $g2, $g3);
if ($n == 3) {
    $gver = 10000 * $g1 + 100 * $g2 + $g3;
} else {
    list($g1, $g2) = sscanf($gretl_version, "%d%c");
    $gver = 10 * $g1 + letter_to_int($g2);
}

$version_data = $pkg . ".versions";
$fp = fopen($version_data, "r") or exit("Unable to open " . $version_data . "!");

$last_version = "";
$got_version = "";

while (!feof($fp)) {
    $pkg_version = trim(fgets($fp));
    if (strlen($pkg_version) > 0) {
        $new_style = 0;
        $n = sscanf($pkg_version, "%d.%d.%d", $p1, $p2, $p3);
	if ($n == 3) {
	    $pver = 10000 * $p1 + 100 * $p2 + $p3;
	} else {
	    list($p1, $p2) = sscanf($pkg_version, "%d%c");
	    $pver = 10 * $p1 + letter_to_int($p2);
	    $new_style = 1;
	}
	if ($pver > $gver) {
	    break;
	} else if ($new_style) {
	    $last_version = "$p1$p2";
	} else {
	    $last_version = "$p1.$p2.$p3";
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

