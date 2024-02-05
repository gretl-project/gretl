# Tests for the Gretl project
Maintainable and readable test code is crucial to establish a good test coverage which in turn enables implementing new features and performing re-factorings without the fear of breaking something.


# Tests
Various types of tests exist. The most common one may be known as unit-tests. Here is quotation of what unit-tests refer to:

	In computer programming, unit testing is a software testing method by which
	individual units of source code—sets of one or more computer program
	modules together with associated control data, usage procedures, and
	operating procedures—are tested to determine whether they are fit for use.
*Source:* https://en.wikipedia.org/wiki/Unit_testing

Unit-tests ensure that a specific function is working as expected. One function might have multiple tests, to catch corner cases or other branches in the code. Having a battery of unit-tests helps to verify that the building blocks of the software work independently from each other.

An alternative kind of tests can be seen the assertion of some 'fundamental' mathematical 'laws'. We would like to build a collection of scripts assertioning that such 'laws' hold when executed using gretl.


# Structure of testing directory
The tests repository has the following structure:

```
tests
	|- practice_scripts
		|- data
	|- test_scripts
		|- commands
        |- data
		|- functions
        |- fundamentals
```

We distinguish between practice scripts and (unit-)tests. This allows one to call either scripts on demand. The test scripts usually do not require much computational time while practice scripts may run much longer which may be inconvenient during the development and building process.

**Practice scripts** are gretl scripts which are supposed to do things as loading a data set, calling some functions and/ or commands as estimating models, conducting some post-estimation analysis *without* asserting any outcome. These scripts can be seen as implicit integration-kind of tests checking that no error is thrown or gretl does not crash when executing such scripts. Users can put any kind of gretl script here -- for instance scripts which let to crashes of gretl in the past. Re-running such scripts helps to make sure that known bugs do not lead to crashes in newer versions any more.

The sub-directory ```./practice_scripts/data``` includes data sets that can be loaded by practice scripts *only*. As data sets may change over time they are not supposed to be used for actual unit-tests. Unit-tests should only use data  which does not change to make sure results are reproducible (see rules 2 and 3 of the *Basic test principles*).

The directory named ```./test_scripts``` comprises actual **test scripts** for gretl commands, functions and some fundamental mathematical laws, respectively. This distinction allows one to execute a battery of existing tests for each type separately: If one works on gretl's C code affecting only gretl functions, it may be unnecessary also to check whether tests for gretl commands run fine (even though it is recommended to execute all tests available at the end before pushing new code to the git master branch).


# Execution of test scripts
The ```./tests``` directory includes the a shell script named `run_tests.sh`.

Print all options by executing `run_tests.sh --help`.

1) ```run_tests.sh --practice```: Execute scripts in ./tests/practice_scripts

2) ```run_tests.sh --commands```: Execute scripts in ./tests/test_scripts/commands

3) ```run_tests.sh --functions```: Execute scripts in ./tests/test_scripts/functions

4) ```run_tests.sh --fundamentals```: Execute scripts in ./tests/test_scripts/fundamentals

5) ```run_tests.sh --all```: Execute options 1 to 4 and run all scripts

All test scripts source the shell-script ```./tests/helper.sh``` which includes all relevant function definitions.


# Write tests

## Basic test principles
1. Write small and specific tests by heavily using helper functions, parametrized tests, gretl's ```assertion()``` function or the ```assertion``` package, not overusing variables, asserting only what’s relevant and avoiding one test for all corner cases.
2. Write self-contained tests by revealing all relevant parameters, insert data right in the test.
3. Write dumb tests by avoiding the reuse of production code and focusing on comparing output values with hard-coded values.
4. KISS > DRY (Don't repeat yourself)
[//]: <> (5. Test close to production by focusing on testing a complete vertical slide and avoiding in-memory databases.)
5. Store test scripts for each function or command as a separate file named according to the rule: ```run_<FUNCTION_NAME>.inp``` or ```run_<COMMAND_NAME>.inp```.


## Given, When, Then
A test should contain three blocks which are separated by one empty line. Each block of code should be as short as possible.

    Print:         Print statement describing which function/ command is tested, and under which conditions.
    Given (Input): Test preparation like creating data or configure mocks
    When (Action): Call the function or action that you like to test
    Then (Output): Execute assertions to verify the correct output or behavior of the action.

Here is an example of code:
```
function void test_stringify (void)
    print "Start testing stringify()."

    # Given
    series y = seq(1, 3)'
    strings S = defarray("A", "B", "C")

    # When
    scalar expected0 = stringify(y, S)
    string expected1 = y[1]
    string expected2 = y[2]
    string expected3 = y[3]

    # Then
    assert(expected0 == FALSE)
    assert(expected1 == "A")
    assert(expected2 == "B")
    assert(expected3 == "C")
end function

test_stringify()
```

See ```./tests/test_scripts/functions/run_stringify.inp``` for some examples

## Function naming
As a convention, the name of the function should start with "test_" followed by

(i) the name or command one wants to test, and

(ii) a description of what are the conditions given.

For instance in ```./tests/test_scripts/run_stringify.inp``` there is the function  named ```test_stringify_nobs_more_nelem()``` which obviously tests gretl's  ```stringify()``` function given the case when the number of observations (nobs) exceeds some number of elements (nelem). Given that gretl functions names cannot be longer than 31 characters, it might be difficult to find self-explanatory function names.


## Use fixed data instead of randomized data
Tests should rely on fixed data only as it can lead to toggling tests which can be hard to debug. Hence, avoid to make use of randomized data but rather pass fixed data to series, matrices, arrays etc.

## Learn from others
Read what other experienced people have to say on testing. Once reference you may start with is https://phauer.com/2019/modern-best-practices-testing-java/.



# References
- https://phauer.com/2019/modern-best-practices-testing-java/
- https://en.wikipedia.org/wiki/Software_testing
