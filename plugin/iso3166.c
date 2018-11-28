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
};

static const struct ccode isocodes[] = {
    {"Afghanistan", "AF", "AFG"},
    {"Åland Islands", "AX", "ALA"},
    {"Albania", "AL", "ALB"},
    {"Algeria", "DZ", "DZA"},
    {"American Samoa", "AS", "ASM"},
    {"Andorra", "AD", "AND"},
    {"Angola", "AO", "AGO"},
    {"Anguilla", "AI", "AIA"},
    {"Antarctica", "AQ", "ATA"},
    {"Antigua and Barbuda", "AG", "ATG"},
    {"Argentina", "AR", "ARG"},
    {"Armenia", "AM", "ARM"},
    {"Aruba", "AW", "ABW"},
    {"Australia", "AU", "AUS"},
    {"Austria", "AT", "AUT"},
    {"Azerbaijan", "AZ", "AZE"},
    {"Bahamas", "BS", "BHS"},
    {"Bahrain", "BH", "BHR"},
    {"Bangladesh", "BD", "BGD"},
    {"Barbados", "BB", "BRB"},
    {"Belarus", "BY", "BLR"},
    {"Belgium", "BE", "BEL"},
    {"Belize", "BZ", "BLZ"},
    {"Benin", "BJ", "BEN"},
    {"Bermuda", "BM", "BMU"},
    {"Bhutan", "BT", "BTN"},
    {"Bolivia (Plurinational State of)", "BO", "BOL"},
    {"Bonaire, Sint Eustatius and Saba", "BQ", "BES"},
    {"Bosnia and Herzegovina", "BA", "BIH"},
    {"Botswana", "BW", "BWA"},
    {"Bouvet Island", "BV", "BVT"},
    {"Brazil", "BR", "BRA"},
    {"British Indian Ocean Territory", "IO", "IOT"},
    {"Brunei Darussalam", "BN", "BRN"},
    {"Bulgaria", "BG", "BGR"},
    {"Burkina Faso", "BF", "BFA"},
    {"Burundi", "BI", "BDI"},
    {"Cabo Verde", "CV", "CPV"},
    {"Cambodia", "KH", "KHM"},
    {"Cameroon", "CM", "CMR"},
    {"Canada", "CA", "CAN"},
    {"Cayman Islands", "KY", "CYM"},
    {"Central African Republic", "CF", "CAF"},
    {"Chad", "TD", "TCD"},
    {"Chile", "CL", "CHL"},
    {"China", "CN", "CHN"},
    {"Christmas Island", "CX", "CXR"},
    {"Cocos (Keeling) Islands", "CC", "CCK"},
    {"Colombia", "CO", "COL"},
    {"Comoros", "KM", "COM"},
    {"Congo (the Democratic Republic of the)", "CD", "COD"},
    {"Congo", "CG", "COG"},
    {"Cook Islands", "CK", "COK"},
    {"Costa Rica", "CR", "CRI"},
    {"Côte d'Ivoire", "CI", "CIV"},
    {"Croatia", "HR", "HRV"},
    {"Cuba", "CU", "CUB"},
    {"Curaçao", "CW", "CUW"},
    {"Cyprus", "CY", "CYP"},
    {"Czechia", "CZ", "CZE"},
    {"Denmark", "DK", "DNK"},
    {"Djibouti", "DJ", "DJI"},
    {"Dominica", "DM", "DMA"},
    {"Dominican Republic", "DO", "DOM"},
    {"Ecuador", "EC", "ECU"},
    {"Egypt", "EG", "EGY"},
    {"El Salvador", "SV", "SLV"},
    {"Equatorial Guinea", "GQ", "GNQ"},
    {"Eritrea", "ER", "ERI"},
    {"Estonia", "EE", "EST"},
    {"Eswatini", "SZ", "SWZ"},
    {"Ethiopia", "ET", "ETH"},
    {"Falkland Islands [Malvinas]", "FK", "FLK"},
    {"Faroe Islands", "FO", "FRO"},
    {"Fiji", "FJ", "FJI"},
    {"Finland", "FI", "FIN"},
    {"France", "FR", "FRA"},
    {"French Guiana", "GF", "GUF"},
    {"French Polynesia", "PF", "PYF"},
    {"French Southern Territories", "TF", "ATF"},
    {"Gabon", "GA", "GAB"},
    {"Gambia", "GM", "GMB"},
    {"Georgia", "GE", "GEO"},
    {"Germany", "DE", "DEU"},
    {"Ghana", "GH", "GHA"},
    {"Gibraltar", "GI", "GIB"},
    {"Greece", "GR", "GRC"},
    {"Greenland", "GL", "GRL"},
    {"Grenada", "GD", "GRD"},
    {"Guadeloupe", "GP", "GLP"},
    {"Guam", "GU", "GUM"},
    {"Guatemala", "GT", "GTM"},
    {"Guernsey", "GG", "GGY"},
    {"Guinea", "GN", "GIN"},
    {"Guinea-Bissau", "GW", "GNB"},
    {"Guyana", "GY", "GUY"},
    {"Haiti", "HT", "HTI"},
    {"Heard Island and McDonald Islands", "HM", "HMD"},
    {"Holy See", "VA", "VAT"},
    {"Honduras", "HN", "HND"},
    {"Hong Kong", "HK", "HKG"},
    {"Hungary", "HU", "HUN"},
    {"Iceland", "IS", "ISL"},
    {"India", "IN", "IND"},
    {"Indonesia", "ID", "IDN"},
    {"Iran (Islamic Republic of)", "IR", "IRN"},
    {"Iraq", "IQ", "IRQ"},
    {"Ireland", "IE", "IRL"},
    {"Isle of Man", "IM", "IMN"},
    {"Israel", "IL", "ISR"},
    {"Italy", "IT", "ITA"},
    {"Jamaica", "JM", "JAM"},
    {"Japan", "JP", "JPN"},
    {"Jersey", "JE", "JEY"},
    {"Jordan", "JO", "JOR"},
    {"Kazakhstan", "KZ", "KAZ"},
    {"Kenya", "KE", "KEN"},
    {"Kiribati", "KI", "KIR"},
    {"Korea (the Democratic People's Republic of)", "KP", "PRK"},
    {"Korea (the Republic of)", "KR", "KOR"},
    {"Kuwait", "KW", "KWT"},
    {"Kyrgyzstan", "KG", "KGZ"},
    {"Lao People's Democratic Republic", "LA", "LAO"},
    {"Latvia", "LV", "LVA"},
    {"Lebanon", "LB", "LBN"},
    {"Lesotho", "LS", "LSO"},
    {"Liberia", "LR", "LBR"},
    {"Libya", "LY", "LBY"},
    {"Liechtenstein", "LI", "LIE"},
    {"Lithuania", "LT", "LTU"},
    {"Luxembourg", "LU", "LUX"},
    {"Macao", "MO", "MAC"},
    {"Macedonia (the former Yugoslav Republic of)", "MK", "MKD"},
    {"Madagascar", "MG", "MDG"},
    {"Malawi", "MW", "MWI"},
    {"Malaysia", "MY", "MYS"},
    {"Maldives", "MV", "MDV"},
    {"Mali", "ML", "MLI"},
    {"Malta", "MT", "MLT"},
    {"Marshall Islands", "MH", "MHL"},
    {"Martinique", "MQ", "MTQ"},
    {"Mauritania", "MR", "MRT"},
    {"Mauritius", "MU", "MUS"},
    {"Mayotte", "YT", "MYT"},
    {"Mexico", "MX", "MEX"},
    {"Micronesia (Federated States of)", "FM", "FSM"},
    {"Moldova (the Republic of)", "MD", "MDA"},
    {"Monaco", "MC", "MCO"},
    {"Mongolia", "MN", "MNG"},
    {"Montenegro", "ME", "MNE"},
    {"Montserrat", "MS", "MSR"},
    {"Morocco", "MA", "MAR"},
    {"Mozambique", "MZ", "MOZ"},
    {"Myanmar", "MM", "MMR"},
    {"Namibia", "NA", "NAM"},
    {"Nauru", "NR", "NRU"},
    {"Nepal", "NP", "NPL"},
    {"Netherlands", "NL", "NLD"},
    {"New Caledonia", "NC", "NCL"},
    {"New Zealand", "NZ", "NZL"},
    {"Nicaragua", "NI", "NIC"},
    {"Niger", "NE", "NER"},
    {"Nigeria", "NG", "NGA"},
    {"Niue", "NU", "NIU"},
    {"Norfolk Island", "NF", "NFK"},
    {"Northern Mariana Islands", "MP", "MNP"},
    {"Norway", "NO", "NOR"},
    {"Oman", "OM", "OMN"},
    {"Pakistan", "PK", "PAK"},
    {"Palau", "PW", "PLW"},
    {"Palestine, State of", "PS", "PSE"},
    {"Panama", "PA", "PAN"},
    {"Papua New Guinea", "PG", "PNG"},
    {"Paraguay", "PY", "PRY"},
    {"Peru", "PE", "PER"},
    {"Philippines", "PH", "PHL"},
    {"Pitcairn", "PN", "PCN"},
    {"Poland", "PL", "POL"},
    {"Portugal", "PT", "PRT"},
    {"Puerto Rico", "PR", "PRI"},
    {"Qatar", "QA", "QAT"},
    {"Réunion", "RE", "REU"},
    {"Romania", "RO", "ROU"},
    {"Russian Federation", "RU", "RUS"},
    {"Rwanda", "RW", "RWA"},
    {"Saint Barthélemy", "BL", "BLM"},
    {"Saint Helena, Ascension and Tristan da Cunha", "SH", "SHN"},
    {"Saint Kitts and Nevis", "KN", "KNA"},
    {"Saint Lucia", "LC", "LCA"},
    {"Saint Martin (French part)", "MF", "MAF"},
    {"Saint Pierre and Miquelon", "PM", "SPM"},
    {"Saint Vincent and the Grenadines", "VC", "VCT"},
    {"Samoa", "WS", "WSM"},
    {"San Marino", "SM", "SMR"},
    {"Sao Tome and Principe", "ST", "STP"},
    {"Saudi Arabia", "SA", "SAU"},
    {"Senegal", "SN", "SEN"},
    {"Serbia", "RS", "SRB"},
    {"Seychelles", "SC", "SYC"},
    {"Sierra Leone", "SL", "SLE"},
    {"Singapore", "SG", "SGP"},
    {"Sint Maarten (Dutch part)", "SX", "SXM"},
    {"Slovakia", "SK", "SVK"},
    {"Slovenia", "SI", "SVN"},
    {"Solomon Islands", "SB", "SLB"},
    {"Somalia", "SO", "SOM"},
    {"South Africa", "ZA", "ZAF"},
    {"South Georgia and the South Sandwich Islands", "GS", "SGS"},
    {"South Sudan", "SS", "SSD"},
    {"Spain", "ES", "ESP"},
    {"Sri Lanka", "LK", "LKA"},
    {"Sudan", "SD", "SDN"},
    {"Suriname", "SR", "SUR"},
    {"Svalbard and Jan Mayen", "SJ", "SJM"},
    {"Sweden", "SE", "SWE"},
    {"Switzerland", "CH", "CHE"},
    {"Syrian Arab Republic", "SY", "SYR"},
    {"Taiwan (Province of China)", "TW", "TWN"},
    {"Tajikistan", "TJ", "TJK"},
    {"Tanzania, United Republic of", "TZ", "TZA"},
    {"Thailand", "TH", "THA"},
    {"Timor-Leste", "TL", "TLS"},
    {"Togo", "TG", "TGO"},
    {"Tokelau", "TK", "TKL"},
    {"Tonga", "TO", "TON"},
    {"Trinidad and Tobago", "TT", "TTO"},
    {"Tunisia", "TN", "TUN"},
    {"Turkey", "TR", "TUR"},
    {"Turkmenistan", "TM", "TKM"},
    {"Turks and Caicos Islands", "TC", "TCA"},
    {"Tuvalu", "TV", "TUV"},
    {"Uganda", "UG", "UGA"},
    {"Ukraine", "UA", "UKR"},
    {"United Arab Emirates", "AE", "ARE"},
    {"United Kingdom of Great Britain and Northern Ireland", "GB", "GBR"},
    {"United States Minor Outlying Islands", "UM", "UMI"},
    {"United States of America", "US", "USA"},
    {"Uruguay", "UY", "URY"},
    {"Uzbekistan", "UZ", "UZB"},
    {"Vanuatu", "VU", "VUT"},
    {"Venezuela (Bolivarian Republic of)", "VE", "VEN"},
    {"Viet Nam", "VN", "VNM"},
    {"Virgin Islands (British)", "VG", "VGB"},
    {"Virgin Islands (U.S.)", "VI", "VIR"},
    {"Wallis and Futuna", "WF", "WLF"},
    {"Western Sahara*", "EH", "ESH"},
    {"Yemen", "YE", "YEM"},
    {"Zambia", "ZM", "ZMB"},
    {"Zimbabwe", "ZW", "ZWE"},
    {NULL, "0", "0"}
};

static const struct ccode fixups[] = {
    {"Vietnam", "VN", "VNM"},
    {"Palestinian Territories", "PS", "PSE"},
    {"Kosovo", "XK", "XKX"}, /* as of 2018, not an ISO country */
    {"South Korea", "KR", "KOR"},
    {"North Korea", "KP", "PRK"},
    {"Great Britain", "GB", "GBR"},
    {"USA", "US", "USA"},
    {NULL, "0", "0"}
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

enum { CNAME = 1, AC2, AC3 };

/* lookup of country name, alpha-2 code or alpha-2 code,
   given one of these and an output specification
*/

char *iso_country (const char *src, int output, int *err)
{
    char *ret = NULL;
    int input = CNAME;
    const struct ccode *ccodes;
    int try_fixups = 1;
    int i, n, match;

    /* what output is requested? */
    if (output != CNAME && output != AC2 && output != AC3) {
	*err = E_INVARG;
	return NULL;
    }

    if (src == NULL || *src == '\0') {
	*err = E_INVARG;
	return NULL;
    }

    /* the string to look up */
    n = strlen(src);
    if (n == 2 && all_upper(src)) {
	input = AC2;
    } else if (n == 3 && all_upper(src)) {
	input = AC3;
    }

    ccodes = isocodes;

 retry:

    for (i=0; ccodes[i].cname != NULL; i++) {
	if (input == CNAME) {
	    match = strncmp(src, ccodes[i].cname, n) == 0;
	} else if (input == AC2) {
	    match = strcmp(src, ccodes[i].ac2) == 0;
	} else {
	    match = strcmp(src, ccodes[i].ac3) == 0;
	}
	if (match) {
	    if (output == CNAME) {
		ret = gretl_strdup(ccodes[i].cname);
	    } else if (output == AC2) {
		ret = gretl_strdup(ccodes[i].ac2);
	    } else {
		ret = gretl_strdup(ccodes[i].ac3);
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
	fprintf(stderr, "iso3166: '%s' was not matched\n", src);
    }

    return ret;
}

gretl_array *iso_country_array (gretl_array *srcs,
				int output,
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
	    result = iso_country(src, output, err);
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
