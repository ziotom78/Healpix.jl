# How to prepare a new release

-   Checkout `master` and pull any changes from the remote repository

-   Update the version number in `CHANGELOG.md` and `Project.toml`

-   Create a commit with comment “Update version to X.Y.Z”

-   Open the commit in the GitHub page and put the following text as a comment:

    ```
    @JuliaRegistrator register

    Release notes:
    BLAH BLAH BLAH
    ```

