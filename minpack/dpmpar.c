double dpmpar_ (int k)
{
    const double ieee_vals[] = {
	2.22044604926e-16,
	2.22507385852e-308,
	1.79769313485e308
    };

    if (k >= 1 && k <= 3) {
	return ieee_vals[k-1];
    } else {
	return 0;
    }
}


