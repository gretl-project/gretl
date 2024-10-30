#include <iostream>
#include <gretl/libgretl.h>
#include <gretl/gretl_win32.h>

#if (defined __GNUC__) && !(defined __clang__)
# define COMPILER_IDENT "GCC " __VERSION__
#else
# define COMPILER_IDENT __VERSION__
#endif

int main ()
{
	gretl_bundle *b = nullptr;
    PRN *prn = nullptr;
	int err = 0;

	libgretl_init();
#if defined(_WIN32)
    win32_cli_read_rc();
#else
    cli_read_rc();
#endif

	std::cout << "compiler used: " << COMPILER_IDENT << std::endl;
	std::cout << "libgretl_version: " << libgretl_version() << std::endl << std::endl;

	b = get_sysinfo_bundle(&err);
	if (!err) {
		prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
		err = gretl_bundle_print(b, prn);
		gretl_print_destroy(prn);
	}

	libgretl_cleanup();
}
