function [NRMSE] = calculateNRMSE(x, y, y_fit)

err = sum((y - y_fit).^2);
RMSE = sqrt(err/length(y));
NRMSE = RMSE/(max(y) - min(y));

end