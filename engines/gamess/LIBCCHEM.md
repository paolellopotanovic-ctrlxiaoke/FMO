## LibCChem Coding Policy

Like GAMESS Fortran, LibCChem seeks to maintain portability and stability. Therefore, the following rules must be followed when submitting new code to LibCChem. Code that does not follow these rules will not be permitted in LibCChem. There will be no exceptions. Note that legacy LibCChem code may not follow these rules; fixing this is a work in progress.

Golden Rule: Make sure that your code compiles with no errors **and no warnings**. The following GCC flags are used for warnings in LibCChem:

```
-\-Wall
-\-Wextra
-\-Wpedantic
-\-Wshadow
-\-Wnon-virtual-dtor
```

Warnings coming from `--Wpedantic` are especially important as these warnings indicate code that is not portable across compilers. Checking warnings across multiple compilers is also highly useful.

Silver Rule: All code submissions and revisions tracking will be done using the GAMESS repository on GitHub. Information on how to submit new features can be found here: [https://github.com/gms-bbg/gamess/wiki](https://github.com/gms-bbg/gamess/wiki)

Rule 1. For stylistic concerns, follow the rules given by the style guide [https://gist.github.com/lefticus/10191322](https://gist.github.com/lefticus/10191322). On top of the rules presented in the style guide, the following extra style rules will also be used:

- `snake_case` is to be used for function definitions and variables, while CamelCase is to be used for class and struct definitions.
- Lines should extend no further than the 80th column to maintain readability.
- The K&R style of bracing is to be used for curly-brace indentation. 

Rule 2. Make sure your code is well commented for readability purposes. Along with the rules in the style guide above, the following comment constructs are to be used to divide your code into subsections:

```
>
>//----------------------//
>//  FIRST-LEVEL HEADER  //
>//----------------------//
>
>//--second-level-header--//
>
```

Rule 3. Code should be written in such a way as to be usable with the GCC `-ffast-math` flag without runtime errors. This enables LibCChem to achieve optimal performance using the `-Ofast` optimization flag.

Rule 4. Whenever possible, use constructs from C++14 instead of constructs from Boost for writing code. Only use constructs from Boost if it does not have a corresponding construct from C++14. Likewise, use constructs from C++14 instead of C-instrinsic constructs whenever possible. For example, use std::arrays instead of C arrays, and use C++14 smart pointers rather than raw pointers. 

Rule 5. Avoid manual memory allocation if possible. If manual memory allocation is required, make sure that every new is matched with a delete operation, and that every malloc is matched with a corresponding free operation.

Rule 6. Make all constructors explicit until implicit construction is required. Additionally, always set the copy and assignment operators equal to delete until they are used.

Rule 7. Whenever a change is made, this change should be recorded in the Git change log included with the submitted code.

Rule 8. While not necessary, it is beneficial for you to scan your code with the cppcheck tool, which can be downloaded [http://cppcheck.sourceforge.net/](http://cppcheck.sourceforge.net/).
