function AnalyseCellCultureData(file_path)
    # Load the CSV file into a DataFrame
    work_dir = pwd();
    #file_path = work_dir * "/Objective_2/Paper2_1/Simulation_Scripts/Cell_Culturing/Data/Lanaro_Browning_Data.csv"
    data = CSV.read(work_dir * file_path, DataFrame);


    # Processing Data
    grouped_data = groupby(data, :PoreSize);

    data_PoreSize_300 = grouped_data[1];
    data_PoreSize_400 = grouped_data[2];
    data_PoreSize_500 = grouped_data[3];
    data_PoreSize_600 = grouped_data[4];

    all_days = [];
    all_mean_density = [];
    all_mean_coverage = [];
    for data_by_PoreSize in groupby(data, :PoreSize);
        days = Float64[];
        mean_density = Float64[];
        mean_coverage = Float64[];
        for data_by_day in groupby(data_by_PoreSize, :Day);
            m_density = sum(data_by_day.S1_Density) / length(data_by_day.S1_Density);
            m_coverage = sum(data_by_day.S2_Coverage) / length(data_by_day.S2_Coverage);
            push!(mean_density, m_density);
            push!(mean_coverage, m_coverage);
            push!(days, data_by_day.Day[1]);
        end
        push!(all_days, days);
        push!(all_mean_density, mean_density);
        push!(all_mean_coverage, mean_coverage);
    end

    avg_density_of_all_exp = sum(hcat(all_mean_density...)[end,:]) / size(hcat(all_mean_density...)[end,:])[1]

    return all_days, all_mean_density, all_mean_coverage, avg_density_of_all_exp

end

function AnalyseBoneFormationVelocity(all_solutions)
    mean_Vn = 0.0
    min_mean_Vn = Inf
    max_mean_Vn = 0.0
    for ii in axes(all_solutions,1)
        sim_mean = mean(mean.(all_solutions[ii].Vₙ))
        if sim_mean < min_mean_Vn
            min_mean_Vn = sim_mean
        end
        if sim_mean > max_mean_Vn
            max_mean_Vn = sim_mean
        end
        mean_Vn += sim_mean
    end
    mean_Vn = mean_Vn / length(all_solutions)
    println("min Vₙ: $min_mean_Vn, avg Vₙ: $mean_Vn, max Vₙ: $max_mean_Vn")
    return (min_mean_Vn, mean_Vn, max_mean_Vn)
end