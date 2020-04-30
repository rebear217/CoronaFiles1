%{
if not(exist('figures','dir'))
    mkdir('figures')
end
if not(exist('synthMatFiles','dir'))
    mkdir('synthMatFiles')    
end
%}

%these do not save PDFs to disk as yet:

%%

%just to show logistic regression under-estimates future data:

close all
OneCountryProjectionAnalysis_LGR(MATdata,'United Kingdom')
%OneCountryProjectionAnalysis_LGR(MATdata,'Italy')
%OneCountryProjectionAnalysis_LGR(MATdata,'Spain')
%OneCountryProjectionAnalysis_LGR(MATdata,'China')
%OneCountryProjectionAnalysis_LGR(MATdata,'Germany')
%OneCountryProjectionAnalysis_LGR(MATdata,'Sweden')

%%

close all
OneCountryProjectionAnalysis_PCA(MATdata,'United Kingdom')

%%

close all
bootstrapODEanalysis(MATdata,'United Kingdom');
%bootstrapODEanalysis(MATdata,'Sweden');
%bootstrapODEanalysis(MATdata,'Italy');
%bootstrapODEanalysis(MATdata,'Spain');
%bootstrapODEanalysis(MATdata,'China');
%bootstrapODEanalysis(MATdata,'Germany');



