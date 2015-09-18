function PlotNestTest( NumTests)
%PLOTNESTTEST Summary of this function goes here
%   Detailed explanation goes here

    data = zeros(NumTests,1);
    errs = zeros(NumTests,1);

    for i = 1:NumTests
        fprintf('%05d: ', i);
        [data(i), errs(i)] = NestTest(i*50+1);
    end

    figure();
    plot((1:NumTests).*50+1, data);
    figure();
    plot((1:NumTests).*50+1, errs);
end
