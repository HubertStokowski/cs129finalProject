%% Calculate fiber to chip coupling for 1550, exactly at 1550.402 nm, average over all powers

% figure();
for ii = 1:length(data_local)
    
    coupling(ii)    = data_local(ii).fiber_to_chip_1550;
    means(ii)       = mean( data_local(ii).Power1550_in_fiber)
    
    wl_loc          = findClosestValue(  data_local(ii).wls, 1549.402  )
    
    vals(ii)        = data_local(ii).Power1550_in_fiber( data_local(ii).wls == wl_loc );
    
%     plot(data_local(ii).wls, data_local(ii).Power1550_in_fiber/max(data_local(ii).Power1550_in_fiber)); hold on;
    
end

%%

% figure(); subplot(311);
% loglog(coupling); hold on;
% subplot(312);
% loglog(means); hold on;
% loglog(power1550); hold on;
% subplot(313);
% loglog(sqrt(0.8*means./power1550)); hold on;
% loglog(sqrt(vals./power1550)); hold on;

coupling_efficiency = mean(sqrt(vals./power1550));

