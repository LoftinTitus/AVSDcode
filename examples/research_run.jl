using AVSDModel

paths = run_research_pipeline()

println("Generated research-grade outputs:")
for key in sort(collect(keys(paths)); by = string)
    println(" - ", key, ": ", paths[key])
end
