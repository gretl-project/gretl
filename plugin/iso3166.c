/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "libgretl.h"
#include "version.h"

struct ccode {
    char *cname; /* Full name of country */
    char ac2[3]; /* Alpha-2 code */
    char ac3[4]; /* Alpha-3 code */
    int nc;      /* Numeric code */
};

/* source: https://www.iso.org/iso-3166-country-codes.html */

static const struct ccode isocodes[] = {
    {"Afghanistan", "AF", "AFG", 4},
    {"Åland Islands", "AX", "ALA", 248},
    {"Albania", "AL", "ALB", 8},
    {"Algeria", "DZ", "DZA", 12},
    {"American Samoa", "AS", "ASM", 16},
    {"Andorra", "AD", "AND", 20},
    {"Angola", "AO", "AGO", 24},
    {"Anguilla", "AI", "AIA", 660},
    {"Antarctica", "AQ", "ATA", 10},
    {"Antigua and Barbuda", "AG", "ATG", 28},
    {"Argentina", "AR", "ARG", 32},
    {"Armenia", "AM", "ARM", 51},
    {"Aruba", "AW", "ABW", 533},
    {"Australia", "AU", "AUS", 36},
    {"Austria", "AT", "AUT", 40},
    {"Azerbaijan", "AZ", "AZE", 31},
    {"Bahamas", "BS", "BHS", 44},
    {"Bahrain", "BH", "BHR", 48},
    {"Bangladesh", "BD", "BGD", 50},
    {"Barbados", "BB", "BRB", 52},
    {"Belarus", "BY", "BLR", 112},
    {"Belgium", "BE", "BEL", 56},
    {"Belize", "BZ", "BLZ", 84},
    {"Benin", "BJ", "BEN", 204},
    {"Bermuda", "BM", "BMU", 60},
    {"Bhutan", "BT", "BTN", 64},
    {"Bolivia (Plurinational State of)", "BO", "BOL", 68},
    {"Bonaire, Sint Eustatius and Saba", "BQ", "BES", 535},
    {"Bosnia and Herzegovina", "BA", "BIH", 70},
    {"Botswana", "BW", "BWA", 72},
    {"Bouvet Island", "BV", "BVT", 74},
    {"Brazil", "BR", "BRA", 76},
    {"British Indian Ocean Territory", "IO", "IOT", 86},
    {"Brunei Darussalam", "BN", "BRN", 96},
    {"Bulgaria", "BG", "BGR", 100},
    {"Burkina Faso", "BF", "BFA", 854},
    {"Burundi", "BI", "BDI", 108},
    {"Cabo Verde", "CV", "CPV", 132},
    {"Cambodia", "KH", "KHM", 116},
    {"Cameroon", "CM", "CMR", 120},
    {"Canada", "CA", "CAN", 124},
    {"Cayman Islands", "KY", "CYM", 136},
    {"Central African Republic", "CF", "CAF", 140},
    {"Chad", "TD", "TCD", 148},
    {"Chile", "CL", "CHL", 152},
    {"China", "CN", "CHN", 156},
    {"Christmas Island", "CX", "CXR", 162},
    {"Cocos (Keeling) Islands", "CC", "CCK", 166},
    {"Colombia", "CO", "COL", 170},
    {"Comoros", "KM", "COM", 174},
    {"Congo (the Democratic Republic of the)", "CD", "COD", 180},
    {"Congo", "CG", "COG", 178},
    {"Cook Islands", "CK", "COK", 184},
    {"Costa Rica", "CR", "CRI", 188},
    {"Côte d'Ivoire", "CI", "CIV", 384},
    {"Croatia", "HR", "HRV", 191},
    {"Cuba", "CU", "CUB", 192},
    {"Curaçao", "CW", "CUW", 531},
    {"Cyprus", "CY", "CYP", 196},
    {"Czechia", "CZ", "CZE", 203},
    {"Denmark", "DK", "DNK", 208},
    {"Djibouti", "DJ", "DJI", 262},
    {"Dominica", "DM", "DMA", 212},
    {"Dominican Republic", "DO", "DOM", 214},
    {"Ecuador", "EC", "ECU", 218},
    {"Egypt", "EG", "EGY", 818},
    {"El Salvador", "SV", "SLV", 222},
    {"Equatorial Guinea", "GQ", "GNQ", 226},
    {"Eritrea", "ER", "ERI", 232},
    {"Estonia", "EE", "EST", 233},
    {"Eswatini", "SZ", "SWZ", 748},
    {"Ethiopia", "ET", "ETH", 231},
    {"Falkland Islands [Malvinas]", "FK", "FLK", 238},
    {"Faroe Islands", "FO", "FRO", 234},
    {"Fiji", "FJ", "FJI", 242},
    {"Finland", "FI", "FIN", 246},
    {"France", "FR", "FRA", 250},
    {"French Guiana", "GF", "GUF", 254},
    {"French Polynesia", "PF", "PYF", 258},
    {"French Southern Territories", "TF", "ATF", 260},
    {"Gabon", "GA", "GAB", 266},
    {"Gambia", "GM", "GMB", 270},
    {"Georgia", "GE", "GEO", 268},
    {"Germany", "DE", "DEU", 276},
    {"Ghana", "GH", "GHA", 288},
    {"Gibraltar", "GI", "GIB", 292},
    {"Greece", "GR", "GRC", 300},
    {"Greenland", "GL", "GRL", 304},
    {"Grenada", "GD", "GRD", 308},
    {"Guadeloupe", "GP", "GLP", 312},
    {"Guam", "GU", "GUM", 316},
    {"Guatemala", "GT", "GTM", 320},
    {"Guernsey", "GG", "GGY", 831},
    {"Guinea", "GN", "GIN", 324},
    {"Guinea-Bissau", "GW", "GNB", 624},
    {"Guyana", "GY", "GUY", 328},
    {"Haiti", "HT", "HTI", 332},
    {"Heard Island and McDonald Islands", "HM", "HMD", 334},
    {"Holy See", "VA", "VAT", 336},
    {"Honduras", "HN", "HND", 340},
    {"Hong Kong", "HK", "HKG", 344},
    {"Hungary", "HU", "HUN", 348},
    {"Iceland", "IS", "ISL", 352},
    {"India", "IN", "IND", 356},
    {"Indonesia", "ID", "IDN", 360},
    {"Iran (Islamic Republic of)", "IR", "IRN", 364},
    {"Iraq", "IQ", "IRQ", 368},
    {"Ireland", "IE", "IRL", 372},
    {"Isle of Man", "IM", "IMN", 833},
    {"Israel", "IL", "ISR", 376},
    {"Italy", "IT", "ITA", 380},
    {"Jamaica", "JM", "JAM", 388},
    {"Japan", "JP", "JPN", 392},
    {"Jersey", "JE", "JEY", 832},
    {"Jordan", "JO", "JOR", 400},
    {"Kazakhstan", "KZ", "KAZ", 398},
    {"Kenya", "KE", "KEN", 404},
    {"Kiribati", "KI", "KIR", 296},
    {"Korea (the Democratic People's Republic of)", "KP", "PRK", 408},
    {"Korea (the Republic of)", "KR", "KOR", 410},
    {"Kuwait", "KW", "KWT", 414},
    {"Kyrgyzstan", "KG", "KGZ", 417},
    {"Lao People's Democratic Republic", "LA", "LAO", 418},
    {"Latvia", "LV", "LVA", 428},
    {"Lebanon", "LB", "LBN", 422},
    {"Lesotho", "LS", "LSO", 426},
    {"Liberia", "LR", "LBR", 430},
    {"Libya", "LY", "LBY", 434},
    {"Liechtenstein", "LI", "LIE", 438},
    {"Lithuania", "LT", "LTU", 440},
    {"Luxembourg", "LU", "LUX", 442},
    {"Macao", "MO", "MAC", 446},
    {"Macedonia (the former Yugoslav Republic of)", "MK", "MKD", 807},
    {"Madagascar", "MG", "MDG", 450},
    {"Malawi", "MW", "MWI", 454},
    {"Malaysia", "MY", "MYS", 458},
    {"Maldives", "MV", "MDV", 462},
    {"Mali", "ML", "MLI", 466},
    {"Malta", "MT", "MLT", 470},
    {"Marshall Islands", "MH", "MHL", 584},
    {"Martinique", "MQ", "MTQ", 474},
    {"Mauritania", "MR", "MRT", 478},
    {"Mauritius", "MU", "MUS", 480},
    {"Mayotte", "YT", "MYT", 175},
    {"Mexico", "MX", "MEX", 484},
    {"Micronesia (Federated States of)", "FM", "FSM", 583},
    {"Moldova (the Republic of)", "MD", "MDA", 498},
    {"Monaco", "MC", "MCO", 492},
    {"Mongolia", "MN", "MNG", 496},
    {"Montenegro", "ME", "MNE", 499},
    {"Montserrat", "MS", "MSR", 500},
    {"Morocco", "MA", "MAR", 504},
    {"Mozambique", "MZ", "MOZ", 508},
    {"Myanmar", "MM", "MMR", 104},
    {"Namibia", "NA", "NAM", 516},
    {"Nauru", "NR", "NRU", 520},
    {"Nepal", "NP", "NPL", 524},
    {"Netherlands", "NL", "NLD", 528},
    {"New Caledonia", "NC", "NCL", 540},
    {"New Zealand", "NZ", "NZL", 554},
    {"Nicaragua", "NI", "NIC", 558},
    {"Niger", "NE", "NER", 562},
    {"Nigeria", "NG", "NGA", 566},
    {"Niue", "NU", "NIU", 570},
    {"Norfolk Island", "NF", "NFK", 574},
    {"Northern Mariana Islands", "MP", "MNP", 580},
    {"Norway", "NO", "NOR", 578},
    {"Oman", "OM", "OMN", 512},
    {"Pakistan", "PK", "PAK", 586},
    {"Palau", "PW", "PLW", 585},
    {"Palestine, State of", "PS", "PSE", 275},
    {"Panama", "PA", "PAN", 591},
    {"Papua New Guinea", "PG", "PNG", 598},
    {"Paraguay", "PY", "PRY", 600},
    {"Peru", "PE", "PER", 604},
    {"Philippines", "PH", "PHL", 608},
    {"Pitcairn", "PN", "PCN", 612},
    {"Poland", "PL", "POL", 616},
    {"Portugal", "PT", "PRT", 620},
    {"Puerto Rico", "PR", "PRI", 630},
    {"Qatar", "QA", "QAT", 634},
    {"Réunion", "RE", "REU", 638},
    {"Romania", "RO", "ROU", 642},
    {"Russian Federation", "RU", "RUS", 643},
    {"Rwanda", "RW", "RWA", 646},
    {"Saint Barthélemy", "BL", "BLM", 652},
    {"Saint Helena, Ascension and Tristan da Cunha", "SH", "SHN", 654},
    {"Saint Kitts and Nevis", "KN", "KNA", 659},
    {"Saint Lucia", "LC", "LCA", 662},
    {"Saint Martin (French part)", "MF", "MAF", 663},
    {"Saint Pierre and Miquelon", "PM", "SPM", 666},
    {"Saint Vincent and the Grenadines", "VC", "VCT", 670},
    {"Samoa", "WS", "WSM", 882},
    {"San Marino", "SM", "SMR", 674},
    {"Sao Tome and Principe", "ST", "STP", 678},
    {"Saudi Arabia", "SA", "SAU", 682},
    {"Senegal", "SN", "SEN", 686},
    {"Serbia", "RS", "SRB", 688},
    {"Seychelles", "SC", "SYC", 690},
    {"Sierra Leone", "SL", "SLE", 694},
    {"Singapore", "SG", "SGP", 702},
    {"Sint Maarten (Dutch part)", "SX", "SXM", 534},
    {"Slovakia", "SK", "SVK", 703},
    {"Slovenia", "SI", "SVN", 705},
    {"Solomon Islands", "SB", "SLB", 90},
    {"Somalia", "SO", "SOM", 706},
    {"South Africa", "ZA", "ZAF", 710},
    {"South Georgia and the South Sandwich Islands", "GS", "SGS", 239},
    {"South Sudan", "SS", "SSD", 728},
    {"Spain", "ES", "ESP", 724},
    {"Sri Lanka", "LK", "LKA", 144},
    {"Sudan", "SD", "SDN", 729},
    {"Suriname", "SR", "SUR", 740},
    {"Svalbard and Jan Mayen", "SJ", "SJM", 744},
    {"Sweden", "SE", "SWE", 752},
    {"Switzerland", "CH", "CHE", 756},
    {"Syrian Arab Republic", "SY", "SYR", 760},
    {"Taiwan (Province of China)", "TW", "TWN", 158},
    {"Tajikistan", "TJ", "TJK", 762},
    {"Tanzania, United Republic of", "TZ", "TZA", 834},
    {"Thailand", "TH", "THA", 764},
    {"Timor-Leste", "TL", "TLS", 626},
    {"Togo", "TG", "TGO", 768},
    {"Tokelau", "TK", "TKL", 772},
    {"Tonga", "TO", "TON", 776},
    {"Trinidad and Tobago", "TT", "TTO", 780},
    {"Tunisia", "TN", "TUN", 788},
    {"Turkey", "TR", "TUR", 792},
    {"Turkmenistan", "TM", "TKM", 795},
    {"Turks and Caicos Islands", "TC", "TCA", 796},
    {"Tuvalu", "TV", "TUV", 798},
    {"Uganda", "UG", "UGA", 800},
    {"Ukraine", "UA", "UKR", 804},
    {"United Arab Emirates", "AE", "ARE", 784},
    {"United Kingdom of Great Britain and Northern Ireland", "GB", "GBR", 826},
    {"United States Minor Outlying Islands", "UM", "UMI", 581},
    {"United States of America", "US", "USA", 840},
    {"Uruguay", "UY", "URY", 858},
    {"Uzbekistan", "UZ", "UZB", 860},
    {"Vanuatu", "VU", "VUT", 548},
    {"Venezuela (Bolivarian Republic of)", "VE", "VEN", 862},
    {"Viet Nam", "VN", "VNM", 704},
    {"Virgin Islands (British)", "VG", "VGB", 92},
    {"Virgin Islands (U.S.)", "VI", "VIR", 850},
    {"Wallis and Futuna", "WF", "WLF", 876},
    {"Western Sahara*", "EH", "ESH", 732},
    {"Yemen", "YE", "YEM", 887},
    {"Zambia", "ZM", "ZMB", 894},
    {"Zimbabwe", "ZW", "ZWE", 716},
    {NULL, "0", "0", 0}
};

static const struct ccode fixups[] = {
    {"Vietnam", "VN", "VNM", 704},
    {"Palestinian Territories", "PS", "PSE", 275},
    {"Kosovo", "XK", "XKX", 383}, /* as of 2018, not an ISO country */
    {"South Korea", "KR", "KOR", 410},
    {"North Korea", "KP", "PRK", 408},
    {"Great Britain", "GB", "GBR", 826},
    {"USA", "US", "USA", 840},
    {"Czech Republic", "CZ", "CZE", 203},
    {"Taiwan Province of China", "TW", "TWN", 158},
    {"Laos", "LA", "LAO", 418},
    {"Ivory Coast", "CI", "CIV", 384},
    {"Congo (Kinshasa)", "CD", "COD", 180},
    {"Congo (Brazzaville)", "CG", "COG", 178},
    {"Hong Kong S.A.R. of China", "HK", "HKG", 344},
    {"Swaziland", "SZ", "SWZ", 748},
    {NULL, "0", "0", 0}
};

static int all_upper (const char *s)
{
    while (*s) {
	if (!isupper(*s)) {
	    return 0;
	}
	s++;
    }

    return 1;
}

enum { CNAME = 1, AC2, AC3, NCD };

/* lookup of country name, alpha-2 code, alpha-3 code
   or numeric code, given one of these and an output
   specification
*/

char *iso_country (const char *src, int output,
		   PRN *prn, int *err)
{
    char *ret = NULL;
    int input = CNAME;
    const struct ccode *ccodes;
    char numout[4];
    int try_fixups = 1;
    int numsrc = 0;
    int i, n, match;

    /* what output is requested? */
    if (output != 0 && output != CNAME &&
	output != AC2 && output != AC3 &&
	output != NCD) {
	*err = E_INVARG;
	return NULL;
    }

    if (src == NULL || *src == '\0') {
	*err = E_INVARG;
	return NULL;
    }

    /* the string to look up */
    n = strlen(src);
    if (n == 3 && isdigit(*src)) {
	input = NCD;
	numsrc = atoi(src);
    } else if (n == 2 && all_upper(src)) {
	input = AC2;
    } else if (n == 3 && all_upper(src)) {
	input = AC3;
    }

    /* automatic output choice? */
    if (output == 0) {
	if (input == CNAME) {
	    output = AC2;
	} else {
	    output = CNAME;
	}
    }

    ccodes = isocodes;

 retry:

    for (i=0; ccodes[i].cname != NULL; i++) {
	if (input == CNAME) {
	    match = strncmp(src, ccodes[i].cname, n) == 0;
	} else if (input == AC2) {
	    match = strcmp(src, ccodes[i].ac2) == 0;
	} else if (input == AC3) {
	    match = strcmp(src, ccodes[i].ac3) == 0;
	} else {
	    match = (numsrc == ccodes[i].nc);
	}
	if (match) {
	    if (output == CNAME) {
		ret = gretl_strdup(ccodes[i].cname);
	    } else if (output == AC2) {
		ret = gretl_strdup(ccodes[i].ac2);
	    } else if (output == AC3) {
		ret = gretl_strdup(ccodes[i].ac3);
	    } else {
		sprintf(numout, "%03d", ccodes[i].nc);
		ret = gretl_strdup(numout);
	    }
	    break;
	}
    }

    if (ret == NULL && try_fixups) {
	try_fixups = 0;
	ccodes = fixups;
	goto retry;
    }

    if (ret == NULL) {
	ret = gretl_strdup("");
	if (prn != NULL) {
	    pprintf(prn, "isocountry: '%s' was not matched\n", src);
	} else {
	    fprintf(stderr, "isocountry: '%s' was not matched\n", src);
	}
    }

    return ret;
}

gretl_array *iso_country_array (gretl_array *srcs,
				int output, PRN *prn,
				int *err)
{
    gretl_array *ret = NULL;
    const char *src;
    char *result;
    int i, n;

    if (gretl_array_get_type(srcs) != GRETL_TYPE_STRINGS) {
	*err = E_TYPES;
	return NULL;
    }

    n = gretl_array_get_length(srcs);
    ret = gretl_array_new(GRETL_TYPE_STRINGS, n, err);

    for (i=0; i<n && !*err; i++) {
	src = gretl_array_get_data(srcs, i);
	if (src == NULL) {
	    *err = E_DATA;
	} else {
	    result = iso_country(src, output, prn, err);
	}
	if (!*err) {
	    *err = gretl_array_set_data(ret, i, result);
	}
    }

    if (*err && ret != NULL) {
	gretl_array_destroy(ret);
	ret = NULL;
    }

    return ret;
}

gretl_array *iso_country_series (const double *x, int n,
				 int output, PRN *prn,
				 int *err)
{
    gretl_array *ret = NULL;
    char src[4];
    char *result;
    int i, k;

    ret = gretl_array_new(GRETL_TYPE_STRINGS, n, err);

    for (i=0; i<n && !*err; i++) {
	k = gretl_int_from_double(x[i], err);
	if (*err || k < 0 || k > 999) {
	    *err = E_INVARG;
	}
	if (!*err) {
	    sprintf(src, "%03d", k);
	    result = iso_country(src, output, prn, err);
	}
	if (!*err) {
	    *err = gretl_array_set_data(ret, i, result);
	}
    }

    if (*err && ret != NULL) {
	gretl_array_destroy(ret);
	ret = NULL;
    }

    return ret;
}
