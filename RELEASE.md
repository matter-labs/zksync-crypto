# Release process

If you want to release the packages on crates.io, follow this process:

1. Install `cargo workspaces`: `cargo install cargo-workspaces`
2. Create a new branch to prepare a release.
3. Change versions in the `Cargo.toml`:
  - `version` in `[workspace.package]`
  - `version` in `[workspace.dependencies]` for all the relevant crates.
4. Run `cargo build`. It must succeed.
5. Commit changes.
6. Run `cargo ws publish --dry-run`. Check the output. It might fail, but it might be OK.
  - `error: config value 'http.cainfo' is not set` can be ignored.
  - There might be warnings, this is OK.
  - There might be errors related to the version resolution, e.g. `failed to select a version`
    (in particular, for `zkevm_test_harness`). It's due to a bug in cargo workspaces.
    Check that the packages it complains about actually have the specified version, and if so,
    it's safe to proceed.
7. Create a PR named `crates.io: Release <version>`. Get a review and merge it.
8. From the main branch _after_ you merge it, run `cargo ws publish --publish-as-is --allow-dirty`.
  - The `--publish-as-is` argument skips the versioning step, which you already did before.
  - The `--allow-dirty` argument is required, because `cargo ws` temporarily removes dev-dependencies
    during publishing.
  - Important: if something fails and you have to do changes to the code, it's safe to run the same
    command again. `cargo ws` will skip already published packages.
9. If something goes wrong, see recommendations below.
10. If everything is OK, create a tag: `git tag v<version>`, e.g. `git tag v0.150.4`
11. `git push --tags`
12. Go to the Releases in the GitHUb, and create a release for published version.
