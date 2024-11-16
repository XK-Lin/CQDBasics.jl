using CQDBasics

cwd = pwd()
subfolers = Dict(
    "Alex 156" => ["13"],
    "Alex 165" => ["59"],
    "Alex 396" => ["1"],
    "FS Low" => ["1"],
    "FS High" => ["1"]
)
for dir ∈ filter(isdir, readdir(cwd, join=true))
    categorize_folders(dir, ("θₙ", "Branching Condition", "Average Method", "Initial μₙ [relative to local B_ext]", "Sigmoid Field", "Bₙ Bₑ Ratio"))
    sort_folders(dir)
    get_subfolder_summary_tables(dir)
end
clean_folders(cwd, false, union(values(subfolers)...))