function countryEstimate = BayesPCA(MATdata,countryStr,Mreduce)
    
    if nargin < 3
        Mreduce = 0;
    end
    
    Nc = length(MATdata.country);
    N = length(MATdata.deathData{1}(1,:));
    
	simTimeLength = MATdata.M;
    MATdata.M = N - Mreduce;
    
    for c = 1:Nc
        MATdata.deathData{c} = MATdata.deathData{c}(:,1:end-Mreduce);
        MATdata.ICaseData{c} = MATdata.ICaseData{c}(:,1:end-Mreduce);        
    end
    
    forceIgnoreList = {'Guinea-Bissau'};    
    countryEstimate = NaN;

    parameters = defaulParameters();
    
    N = length(MATdata.country);
    PCAdims = 5;
    minDeathsToIncludeCountry = parameters.limitDeaths;
    movingAverage = parameters.movingAverageDays;
        
    cList = [];
    j = 1;
    
    compactData = {};
    DisplayOff = true;
    
    for ctry = 1:N
        if ~strcmp(MATdata.country{ctry},forceIgnoreList)
            Ddata = sum(MATdata.deathData{ctry},1);
            if not(strcmp('Global',MATdata.country{ctry})) && ...
                (sum(Ddata) > minDeathsToIncludeCountry)
            
                [PCD,y,x] = getPCD(Ddata,MATdata.country{ctry},DisplayOff);
                cList = [cList ctry];

                compactData{j} = y/max(y);
                j = j + 1;
            end
        end
    end    
    
    lenmx = 1;
    for c = 1:length(cList)
        l = length(compactData{c});
        lenmx = max(l,lenmx);
    end
    
    NN = NaN(length(cList),lenmx);
    smallTimes = 1:min([simTimeLength lenmx]);
    
    for c = 1:length(cList)
        NN(c,1:length(compactData{c})) = compactData{c};
    end
    
    [C, ss, M, X,Ye] = ppca_mv(NN,PCAdims,0);
        
    pl = zeros(1,length(cList));
    pl2 = zeros(1,length(cList));
    
	iC = 1;
    
    for c = 1:length(cList)
        ctry = cList(c);
        
        if movingAverage > 0
            Y = movmean(Ye(c,smallTimes),movingAverage);
        else
            Y = Ye(c,smallTimes);
        end

        Ddata = sum(MATdata.deathData{ctry},1);
        [PCD,y,x] = getPCD(Ddata,MATdata.country{ctry});
        Y = max(y)*Y;
        CY = cumsum(Y);
        
        if strcmp(MATdata.country{ctry},countryStr)
            countryEstimate = CY;
        end
        
    end
    
end