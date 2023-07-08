using Documenter, AutoVectors

makedocs(sitename="My Documentation")

deploydocs(
    root = "<current-directory>",
    target = "build",
    dirname = "",
    repo = "<required>",
    branch = "gh-pages",
    deps = nothing | <Function>,
    make = nothing | <Function>,
    devbranch = nothing,
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#", devurl => devurl],
    forcepush = false,
    deploy_config = auto_detect_deploy_system(),
    push_preview = false,
    repo_previews = "github.com/srwhite59/AutoVectors.jl.git",
    branch_previews = "gh-pages",
)
