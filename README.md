FastCA
==============================================================================
This is an effective and efficient tool for combinatorial interactive test. It can generate smaller covering array in a much shorter time, compared to other solvers.


Prerequisites
--------
- Java Runtime Environment
- g++ with C++11 supported

Build
-------
```shell
cd fastca
make
```
The executable file of *FastCA* will appear under the `fastca` directory.

Usage
-----
The command for running *FastCA* is:

```shell
./FastCA --file \<filename1\> --strength <covering strength> --time <time budget> --outfile <filename2>
```

The model file of the software under testing is described in `filename1` (see below).

The covering strength is an integer indicating the testing requirement. It means every combination of any `t` parameter values in the software must be tested at least once.

The time budget is given in CPU seconds. Once the time budget is reached, *FastCA* prints the generated covering array in the output file `filename2`.

Example
------
The model of the software under test is illustrated as the following example.

```
[System]
Name: Example
[Parameter]
P_0(int): 0, 1, 2, 3
P_1(enum): x, y
P_2(int): 0, 1, 2
[Constraint]
P_0=0 || P_1=y 
P_0!=0 || P_1=x || P_2!=1 
[Test Set]
P_0, P_1, P_2
0, y, 2
*, x, 2
```

The name of SUT is specified in the *System* section. The name, type,
and possible values of each parameter are given in the *Parameter*
section. Currently, *FastCA* supports the integer and enumerated types.
Other types should be easily presented as enumerated.

The *Constraint* section specifies the constraints on the permissible
combination of values. It is written as a conjunction of disjunctions
over equality and inequality terms, one disjunction clause per line. The
solver will generate a covering array subjecting to all constraints
indicated here.

*FastCA* also supports constructing covering arrays
extended from given patterns which are specified in the *Test Set*
section. The first line of *Test Set* indicates the names of concerned
parameters, while the following lines describe the patterns which are
required to appear at least once in the generated solution. The symbol
‘\*’ means $any$ value. These patterns are also called *seeds*. The
*Constraint* and *Test Set* sections are both optional.

Once the given cutoff time is reached, *FastCA* prints the generated CA,
of which each line is a complete assignment to the parameters. The
parameters and the set of their values are also printed to indicate the
order of parameters in the generated CA. 

For the example above, *FastCA* will generate a 3-way CA of size 14, and the first configuration is (P0=0, P1=y, P2=2). All valid 3-way interactions shall be covered
by testing the SUT with the generated CA.

```
Parameters:
P_0: [0, 1, 2, 3]
P_1: [x, y]
P_2: [0, 1, 2]
Configurations:
1st   0 y 2
2nd   0 x 1
3rd   0 x 0
…
14th   3 y 2
Found Covering Array of size : 14

```

Possible Extension
------------------

For industrial or safety-critical systems, it usually requires higher covering strength testing for subsets of components, while keeping relatively lower covering strength for other components.
In this case, variable strength covering array may be more effective and efficient.

Currently, *FastCA* only supports generating CAs with uniform covering strength for all options. But it should be easily extended to support variable covering strength by maintaining independent data structures for score computation of each covering strength.

Summary of Experimental results
----
Please refer to the file `tables.pdf`.

- Table 1-3 present the score computation of `TCA` for solving CCAG with covering stenght between 2-4.

- Table 4 and 5 compare `TCA-opt` with `TCA` for 2-way and 3-way CCAG.

- Table 6-8 compare `FastCA` against its competitors for 2-4 way CCAG.

- Table 9 present the experimental results about alternative versions of `FastCA` on 3-way real-world instances.

Additional tables with cutoff time of 5 hours
----
Please refer to the file `additional_table.pdf`
