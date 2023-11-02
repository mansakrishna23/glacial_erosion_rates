%% Modelling time scale biases in erosion rates
% This script attempts to replicate the numerical model described by <https://www.science.org/doi/10.1126/sciadv.1600204 
% Ganti et al. (2016)>, in which they found time scale biases in erosion rates 
% could be generated by random independent time series of hiatuses and erosional 
% pulses generated by a heavy tailed distribution
%% Parameters
% Parameter values as defined by Ganti et al. (2016)

tail_index = 0.5; % pareto parameter
upper_bound = 200e3; % max hiatus length [y]
t_max = 5e6; % max time [y]
erosion_mag = 10; % erosion magnitude [mm]
N = 1000; % number of independent model runs

% cutoff to ensure that start of time series is independent of hiatus length
cutoff = 500e3; 
% time scales at which to evalusate erosion rates
time_scales = sort([10.^(0:6), upper_bound,4e5, 2e6,4e6]);
%% Example model run
% Representative time series of erosion magnitude and time before present

time_bp = sample_pareto(1, tail_index, upper_bound, t_max, cutoff);

figure(1)
clf
hold on
set(gca,'FontSize',12)
time_erosion = [time_bp,time_bp,time_bp]';
sz = size(time_bp);
erosion_pulses = [zeros(sz),erosion_mag*ones(sz),zeros(sz)]';
plot(time_erosion(:), erosion_pulses(:),'color','k','LineWidth',1)

xlim([0 t_max])
ylim([0 12])

xlabel('Time BP [My]')
ylabel('Erosion [mm]')
%% Run model
% This implementation attempts to follow Ganti et al. (2016) as described

% erosion rates for given time scales
erosion_rates = run_model(N, tail_index, upper_bound, t_max, cutoff, erosion_mag, time_scales);
% erosion rates ignoring periods with 0 total erosion
erosion_rates_nonzero = erosion_rates;
erosion_rates_nonzero(erosion_rates==0) = nan;
% counts of zero and non zero erosion rates for given time scales
count_zeros = sum(erosion_rates == 0,1);
count_nonzeros = sum(erosion_rates > 0,1);

% FIGURE: model results
figure(2)
clf

% SUBPLOT: apparent erosion rates by time scale
subplot(3,3,1:6)
hold on
box on
set(gca,'XScale','log','yscale','log','FontSize',10,'layer','top')
xlim([10^-0.5,1e7])
ylim([1e-4,1e3])
xticks(10.^(0:7))
yticks(10.^(-5:2))
xlabel('Time scale [y]')
ylabel('Erosion rate [mm/y]')

% plot upper bound
patch([upper_bound,upper_bound,1e7,1e7],[ylim,flip(ylim)],0.9*[1 1 1],'edgecolor','none')
line([upper_bound,upper_bound],ylim,'color','k','linewidth',0.75,'linestyle','--')

% plot data and means
s = scatter(time_scales,erosion_rates, ...
    10,'r','filled','MarkerFacealpha',0.1);
s_1 = plot(time_scales,mean(erosion_rates,1,'omitmissing'),'ok', ...
    'markersize',7.5,'markerfacecolor','k','DisplayName','All');
s_2 = plot(time_scales,mean(erosion_rates_nonzero,1,'omitmissing'),'or', ...
    'markersize',7.5,'markerfacecolor','r','DisplayName','Non-zero');

% linear fit for mean erosion rates
con = time_scales<upper_bound;
x = time_scales(con);
y = erosion_rates(:,con);
y_m = mean(y,1,'omitmissing');
p = polyfit(log10(x),log10(y_m),1);
f = 10.^polyval(p,log10(x));
plot(x,f,'k','linewidth',1)
text(2*1e2,1e-2,sprintf('%1.2f',p(1)),'Color','k')

% linear fit for non zero means
y = erosion_rates_nonzero(:,con);
y_m = mean(y,1,'omitmissing');
p = polyfit(log10(x),log10(y_m),1);
f = 10.^polyval(p,log10(x));
plot(x,f,'r','linewidth',1)
text(2*1e2,2*1e0,sprintf('%1.2f',p(1)),'Color','r')


% legend
s_l1 = plot(nan,nan,'o','markerfacecolor',[0.5 0.5 0.5], 'markeredgecolor','none', 'markersize',7.5,'DisplayName','Mean');
s_l2 = plot(nan,nan,'o','markerfacecolor',[0.8 0.8 0.8], 'markeredgecolor','none', 'markersize',3,'DisplayName','Data');
legend([s_1,s_2,s_l1,s_l2],'location','northwest','NumColumns',2)

% SUBPLOT: zero and non-zero data count
subplot(3,3,7:9)
hold on
box on
set(gca,'FontSize',10,'layer','top','xAxisLocation','top')
ylim([0 N])
ylabel('Count')
xlim([-0.5 7])
xlabel([])
xticklabels([])

patch([log10(upper_bound),log10(upper_bound),7,7],[ylim,flip(ylim)],0.9*[1 1 1],'edgecolor','none')
b = bar(log10(time_scales),[count_nonzeros;count_zeros],'stacked', 'FaceColor','flat','EdgeColor','none');
b(1).CData = [1 0 0];
b(2).CData = [0 0 0];
line([log10(upper_bound),log10(upper_bound)],ylim,'color','k','linewidth',0.75,'linestyle','--')

%% Functions

function time_bp = sample_pareto(seed, tail_index, upper_bound, t_max, cutoff)
    % pull random samples of hiatuses from truncated pareto distribution
    rng(seed) % For reproducibility
    N = t_max/100; % number of samples
    u = rand(N, 1); % N uniform random numbers (0, 1)
    x_tp = round((u).^(-1/tail_index)); % pareto samples
    x_tp = x_tp(x_tp<upper_bound); % truncate pareto

    % build time series (0, t_max)
    time_bp = cumsum(x_tp) - cutoff;
    ix = time_bp>0 & time_bp<=t_max;
    time_bp = time_bp(ix);
end

function erosion_rates = run_model(N, tail_index, upper_bound, t_max, cutoff, erosion_mag, time_scales)
    erosion_rates = zeros(N,length(time_scales)); % to store erosion rates for each time scale for each model run
    
    for i = 1:N
        time_bp = sample_pareto(i, tail_index, upper_bound, t_max,cutoff);
        erosion = ones(size(time_bp)) * erosion_mag;
        erosion_total = cumsum(erosion);
        ix = sum(time_bp <= time_scales, 1);
        non_zero = ix>0;
        rates = erosion_total(ix(non_zero))./time_bp(ix(non_zero));
        erosion_rates(i, non_zero) = rates;
    end
end
