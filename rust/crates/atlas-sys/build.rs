const ATLAS_VERSION: &str = "0.46.0";

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=ATLAS_DIR");
    println!("cargo:rerun-if-env-changed=DOCS_RS");

    if bindman_utils::is_docs_rs() {
        return;
    }

    bindman_utils::validate_build_mode(cfg!(feature = "system"), cfg!(feature = "vendored"));

    if cfg!(feature = "system") {
        build_system();
    } else {
        build_vendored();
    }
}

#[cfg(feature = "system")]
fn build_system() {
    let (root, include, lib_dir) = bindman_utils::cmake_find_package("atlas", ATLAS_VERSION);

    println!("cargo:rustc-link-search=native={}", lib_dir.display());
    println!("cargo:rustc-link-lib=dylib=atlas");

    println!("cargo:root={}", root.display());
    println!("cargo:include={}", include.display());
}

#[cfg(not(feature = "system"))]
fn build_system() {
    unreachable!("build_system called without system feature");
}

#[cfg(feature = "vendored")]
fn build_vendored() {
    use std::env;
    use std::fs;
    use std::path::PathBuf;
    use std::process::Command;

    const ECBUILD_REPO: &str = "https://github.com/ecmwf/ecbuild.git";
    const ECBUILD_TAG: &str = "3.13.1";
    const ATLAS_REPO: &str = "https://github.com/ecmwf/atlas.git";

    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR not set"));
    let src_dir = out_dir.join("src");
    let build_dir = out_dir.join("build");
    let install_dir = out_dir.join("install");

    fs::create_dir_all(&src_dir).expect("Failed to create src directory");
    fs::create_dir_all(&build_dir).expect("Failed to create build directory");

    let eckit_root =
        env::var("DEP_ECKIT_SYS_ROOT").expect("DEP_ECKIT_SYS_ROOT not set - eckit-sys dependency");

    let ecbuild_src = bindman_utils::git_clone(ECBUILD_REPO, ECBUILD_TAG, &src_dir.join("ecbuild"));
    let atlas_src = bindman_utils::git_clone(ATLAS_REPO, ATLAS_VERSION, &src_dir.join("atlas"));

    let ecbuild_bin = ecbuild_src.join("bin/ecbuild");
    let num_jobs = bindman_utils::build_parallelism();

    let mut cmd = Command::new(&ecbuild_bin);
    cmd.current_dir(&build_dir)
        .arg(format!("--prefix={}", install_dir.display()))
        .arg("--")
        .arg(&atlas_src)
        .arg(format!("-DCMAKE_PREFIX_PATH={eckit_root}"))
        .arg(format!(
            "-DCMAKE_BUILD_TYPE={}",
            bindman_utils::cmake_build_type()
        ))
        .arg("-DENABLE_TESTS=OFF")
        .arg("-DBUILD_TESTING=OFF")
        .arg("-DENABLE_DOCS=OFF")
        .arg("-DENABLE_FORTRAN=OFF")
        .arg("-DENABLE_SANDBOX=OFF")
        .arg("-DENABLE_CLANG_TIDY=OFF")
        .arg(format!(
            "-DENABLE_OMP={}",
            bindman_utils::on_off(cfg!(feature = "omp"))
        ));

    #[cfg(target_os = "macos")]
    cmd.arg("-DCMAKE_INSTALL_NAME_DIR=@rpath");

    bindman_utils::run_command(&mut cmd, "ecbuild configure atlas");

    bindman_utils::run_command(
        Command::new("cmake")
            .args(["--build", ".", "--parallel", &num_jobs])
            .current_dir(&build_dir),
        "cmake build atlas",
    );

    bindman_utils::run_command(
        Command::new("cmake")
            .args(["--install", "."])
            .current_dir(&build_dir),
        "cmake install atlas",
    );

    let lib_dir = bindman_utils::resolve_lib_dir(&install_dir);

    println!("cargo:rustc-link-search=native={}", lib_dir.display());
    println!("cargo:rustc-link-lib=dylib=atlas");
    bindman_utils::link_cpp_stdlib();

    println!("cargo:root={}", install_dir.display());
    println!("cargo:include={}", install_dir.join("include").display());
}

#[cfg(not(feature = "vendored"))]
fn build_vendored() {
    unreachable!("build_vendored called without vendored feature");
}
