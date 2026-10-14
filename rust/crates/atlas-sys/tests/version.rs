//! Guard against drift between the crate version and the repo's `VERSION`
//! file (cargo cannot read it dynamically).

#[test]
fn crate_version_matches_repo_version_file() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../VERSION");
    let Ok(version) = std::fs::read_to_string(path) else {
        // Not building from the repo checkout (e.g. a packaged crate).
        return;
    };
    assert_eq!(
        version.trim(),
        env!("CARGO_PKG_VERSION"),
        "rust/Cargo.toml workspace version and VERSION file are out of sync"
    );
}
