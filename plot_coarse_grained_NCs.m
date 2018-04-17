function plot_coarse_grained_NCs(distn,ts,omit_oocyte)
%created 11/5/16
%last edited 2/10/17
%JH
%plots data for distn of mRNA in nurse cells onto a representation of NCs
%see caceres development paper 2005 for NCs map
%based on previous version of this fn in coarse-graining project folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%distn is a matrix T x 16 containing rna levels for all the cells at different time pts as rows
%ts is a vector of length T containing the corresponding times
%omit_oocyte is a logical giving whether the number of particles in the oocyte should be plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
if length(ts) ~= size(distn,1)
error('input data wrong size');
end

%map integers 1:16 onto their position in NCs layout diagram
NCs_map = [2,3;
    2,2;
    3,3;
    3,2;
    2,4;
    2,1;
    3,4;
    3,1;
    1,3;
    1,2;
    4,3;
    4,2;
    1,4;
    1,1;
    4,4;
    4,1];

connections{1} = [2,3,5,9];
connections{2} = [1,4,6,10];
connections{3} = [1,7,11];
connections{4} = [2,8,12];
connections{5} = [1,13];
connections{6} = [2,14];
connections{7} = [3,15];
connections{8} = [4,16];
connections{9} = 1;
connections{10} = 2;
connections{11} = 3;
connections{12} = 4;
connections{13} = 5;
connections{14} = 6;
connections{15} = 7;
connections{16} = 8;

map_reverted = fliplr(NCs_map);
for j=1:16
    for k=1:length(connections{j})
        z = connections{j}(k);
        coords(j,:) = [map_reverted(j,:),map_reverted(z,:)]-0.5;
    end
end

M = zeros(4); %store mRNA concentration (proportion of total mRNA) in M
for i=1:size(distn,1)
    %loop over the cells
    for j=1:16
        M(NCs_map(j,1),NCs_map(j,2)) = distn(i,j);
    end
    figure(i);
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    if ~isnan(ts)
        title(sprintf('Nurse cells at t=%2.1f ',ts(i)));
    end
    hold on;
    if omit_oocyte;
        M(NCs_map(1,1),NCs_map(1,2))=0; %remove oocyte for plotting
        mm=max(max(M));
    else
        mm=max(max(distn));
    end
    for i1=1:4
        for i2=1:4
            rectangle('Position',[i1-1,i2-1,1,1],...
                'Curvature',[1,1], 'FaceColor',M(5-i2,i1)/mm*ones(1,3)); %plots/counts from bottom left, rather than top left
        end
    end
    for j=1:16
        line(coords(j,[1,3]),coords(j,[2,4]),'Color','r');
    end
%     %for j=1:16
%         A = NCs_map(j,:)
%         B = NCs_map(connections{j},:)
%     [4.5-B(1,1),4.5-A(1,1)]
%     [4.5-B(1,2),4.5-A(1,2)]
%         line([4.5-B(1,1),A(1,1)-0.5],[4.5-B(1,2),A(1,2)-0.5],'Color','r');
%         
%     %end
    
    %     figure;
    %     rectangle('Position',[NCs_map(j,1)-1,NCs_map(j,2)-1,1,1],...
    %   'Curvature',[1,1], 'FaceColor','r')
end
