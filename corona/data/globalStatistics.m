countries = length(MATdata.country);

totalDeaths = zeros(1,MATdata.M);
totalIs = zeros(1,MATdata.M);

globalPopulation = 0;

for c = 1:countries
    thisDeaths = sum(MATdata.deathData{c},1);
    totalDeaths = totalDeaths + thisDeaths;

	thisI = sum(MATdata.ICaseData{c},1);
    totalIs = totalIs + thisI;

    if isnan(MATdata.population{c})
        MATdata.population{c} = [];
    end
    if ~isempty(MATdata.population{c})
        globalPopulation = globalPopulation + MATdata.population{c};
    end
end

MATdata.country{countries+1} = 'Global';
MATdata.deathData{countries+1} = totalDeaths;
MATdata.ICaseData{countries+1} = totalIs;
MATdata.population{countries+1} = globalPopulation;
