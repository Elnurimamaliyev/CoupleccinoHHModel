
%%
data2 = load ('20240701-B-M-1');
%%
% pre_cell_15
% post_cell_15
% pre_cell_pharma_12
% post_cell_pharma_12
data1 = load('20240701-S-A-1');
% pre_pharma(1:15)
data1.pre_cell_15

data3 = load('20240701-B-M-2');



%% Control
figure;

subplot(1,2,1)
plot(data1.pre_cell_15(15, 239001:259001))
% data_to_plot = np.concatenate((pre_cell_11[7, 179000:199001], pre_cell_11[7, 239001:259001]))
% xlim([179000:199000, 239001:259000])
subplot(1,2,2)
plot(data1.post_cell_15(3, 239001:259001))

% Pharma
figure;
subplot(1,2,1)
plot(data1.pre_cell_pharma_28(3, 239001:259001))
subplot(1,2,2)
plot(data1.post_cell_pharma_28(3, 239001:259001))

%% Control
figure;

subplot(1,2,1)
% 19,22,23
plot(data3.pre_cell_pharma_30(22, 239001:259001))
subplot(1,2,2)
plot(data3.post_cell_pharma_30(22, 239001:259001))


batmanreal.resp = data3.pre_cell_pharma_30([19,22,23],:);
% % Pharma
% figure;
% subplot(1,2,1)
% plot(data2.pre_cell_15(4, 239001:259001))
% subplot(1,2,2)
% plot(data2.post_cell_15(4, 239001:259001))

%%
figure;

num_recordings = size(data2.pre_cell_15, 1); % Number of recordings (assuming pre_cell_15 and post_cell_15 have the same number of recordings)

for i = 1:num_recordings
    % Plot for pre_cell_15
    subplot(num_recordings, 2, 2*i - 1);
    plot(data2.pre_cell_15(i, 239001:259001));
    title(['Pre Cell Recording ', num2str(i)]);
    xlabel('Time');
    ylabel('Potential');
    
    % Plot for post_cell_15
    subplot(num_recordings, 2, 2*i);
    plot(data2.post_cell_15(i, 239001:259001));
    title(['Post Cell Recording ', num2str(i)]);
    xlabel('Time');
    ylabel('Potential');
end

%%
figure;

num_recordings = size(data3.pre_cell_15, 1); % Number of recordings (assuming pre_cell_15 and post_cell_15 have the same number of recordings)

for i = 1:num_recordings
    % Plot for pre_cell_15
    subplot(num_recordings, 2, 2*i - 1);
    plot(data3.pre_cell_15(i, 239001:259001));
    title(['Pre Cell Recording ', num2str(i)]);
    xlabel('Time');
    ylabel('Potential');
    
    % Plot for post_cell_15
    subplot(num_recordings, 2, 2*i);
    plot(data3.post_cell_15(i, 239001:259001));
    title(['Post Cell Recording ', num2str(i)]);
    xlabel('Time');
    ylabel('Potential');
end

%%
figure;

num_recordings = size(data3.pre_cell_pharma_30, 1); % Number of recordings (assuming pre_cell_15 and post_cell_15 have the same number of recordings)


for i = 1:num_recordings
    % Plot for pre_cell_15
    subplot(num_recordings, 2, 2*i - 1);
    plot(data3.pre_cell_pharma_30([19,22,23], 239001:259001));
    title(['Pre Cell Recording ', num2str(i)]);
    xlabel('Time');
    ylabel('Potential');
    
    % Plot for post_cell_15
    subplot(num_recordings, 2, 2*i);
    plot(data3.post_cell_pharma_30(i, 239001:259001));
    title(['Post Cell Recording ', num2str(i)]);
    xlabel('Time');
    ylabel('Potential');
end

19,22,23
%%
% Cell 1 - pre pharma

% Cell 15 - after pharma
data2.post_cell_pharma_12(10, 239001:259001)

%%

figure;

num_recordings = size(data2.post_cell_pharma_12, 1); % Number of recordings (assuming pre_cell_15 and post_cell_15 have the same number of recordings)

for i = 1:num_recordings
    % Plot for pre_cell_15
    subplot(num_recordings, 2, 2*i - 1);
    plot(data2.pre_cell_pharma_12(i, 239001:259001));
    title(['Pre Cell Recording ', num2str(i)]);
    xlabel('Time');
    ylabel('Potential');
    
    % Plot for post_cell_15
    subplot(num_recordings, 2, 2*i);
    plot(data2.post_cell_pharma_12(i, 239001:259001));
    title(['Post Cell Recording ', num2str(i)]);
    xlabel('Time');
    ylabel('Potential');
end
