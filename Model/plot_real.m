%% import Data
% Names

Data_big.Data1  = load('20240628-B-M-1');
Data_big.Data2  = load('20240628-B-M-2');
Data_big.Data3  = load('20240628-B-M-3');
Data_big.Data4  = load('20240628-S-A-1');
Data_big.Data5  = load('20240628-S-A-3');
Data_big.Data6  = load('20240701-S-A-1');
Data_big.Data7  = load('20240701-B-M-1');
Data_big.Data8  = load('20240701-B-M-2');
Data_big.Data9  = load('20240702-B-M-1');
Data_big.Data10 = load('20240702-B-M-2');
Data_big.Data11 = load('20240703-B-M-2');
Data_big.Data12 = load('20240703-B-M-4');
Data_big.Data13 = load('20240703-S-A-2');
Data_big.Data14 = load('20240703-S-A-3');
%%
Save(Data_big)
%%
figure(1);
subplot(2,1,1);
plot(Data_big.Data1.pre_cell_13(5,:));
hold on 
plot(Data_big.Data1.pre_cell_pharma_29(20,:));

subplot(2,1,2);
plot(Data_big.Data1.post_cell_13(5,:))
hold on 
plot(Data_big.Data1.post_cell_pharma_29(20,:));

%%
figure(2);
subplot(2,1,1);
plot(Data_big.Data2.pre_cell_15(5,:));
hold on 
plot(Data_big.Data2.pre_cell_pharma_38(20,:));

subplot(2,1,2);
plot(Data_big.Data2.post_cell_15(5,:))
hold on 
plot(Data_big.Data2.post_cell_pharma_38(20,:));

%%
figure(3);
subplot(2,1,1);
plot(Data_big.Data3.pre_cell_15(11,:));
hold on 
plot(Data_big.Data3.pre_cell_pharma_15(10,:));

subplot(2,1,2);
plot(Data_big.Data3.post_cell_15(11,:))
hold on 
plot(Data_big.Data3.post_cell_pharma_15(10,:));
%%
figure(4);
subplot(2,1,1);
plot(Data_big.Data4.pre_cell_12(11,:));
hold on 
plot(Data_big.Data4.pre_cell_pharma_9(5,:));

subplot(2,1,2);
plot(Data_big.Data4.post_cell_12(11,:))
hold on 
plot(Data_big.Data4.post_cell_pharma_9(5,:));

%%
figure(5);
subplot(2,1,1);
plot(Data_big.Data5.pre_cell_15(11,:));
hold on 
plot(Data_big.Data5.pre_cell_pharma_15(10,:));

subplot(2,1,2);
plot(Data_big.Data5.post_cell_15(11,:))
hold on 
plot(Data_big.Data5.post_cell_pharma_15(10,:));

%%
figure(6);
subplot(2,1,1);
plot(Data_big.Data6.pre_cell_15(11,:));
hold on 
plot(Data_big.Data6.pre_cell_pharma_15(10,:));

subplot(2,1,2);
plot(Data_big.Data6.post_cell_15(11,:))
hold on 
plot(Data_big.Data6.post_cell_pharma_15(10,:));
%%

