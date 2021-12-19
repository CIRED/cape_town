%importation des fichiers de temps de transport envoy?s par la ville du Cap
%juste besoin de le lancer une fois pour que ?a cr?e les fichiers .mat,
%puis apr?s on peut faire des simulations sans avoir besoin de relancer ?a
%?a tourne en 3 ou 4 minutes

%liste des fichiers: tt c'est pour travel time et t sh c'ets pour transport
%share, i.e. la matrice OD
liste={'bustsh','carstt','taxitsh','totaltsh','traintt','bustt','cartsh','taxitt','traintsh','walktsh'};

%on importe les fichiers
for i=liste
    i=char(i);
    eval([i,'=importfile_mat(''./data/transport/',i,'.csv'');'])
end

%importation et sauvegarde de la liste des codes des zones
liste1=importfile_mat('./data/transport/carstt.csv',1, 1);
liste2=liste1;
%save('liste.mat','liste1','liste2');

%% importation des coordonn?es des zones

%fichier de donn?es sur els zones: emplois et coordonn?es
[X,Y,TZ2013,Area,diss,Zones,BY5Origins,BY5Destina,PTODOrigin,PTODDestin,emploi_TAZ,emploi_T_1,emploi_T_2,emploi_T_3,emploi_T_4,emploi_T_5,emploi_T_6,emploi_T_7,emploi_T_8,job_total,job_dens,Ink1,Ink2,Ink3,Ink4,job_tot] ...
    = importfile_TAZ('./data/TAZ_amp_2013_proj_centro2.csv');

%on r?-ordonne les temps de transport comme il faut, pour qu'ils soient
%dans le m?me ordre que les zones
cars=zeros(size(TZ2013));
bus=zeros(size(TZ2013));
taxi=zeros(size(TZ2013));
train=zeros(size(TZ2013));
for index1=1:length(TZ2013),
   for index2=1:length(TZ2013),
       choix1=(liste1==TZ2013(index1));
       choix2=(liste2==TZ2013(index2));
       if (sum(choix1)>0)&&(sum(choix2)>0),
           cars(index1,index2)=carstt(choix1,choix2);
           bus(index1,index2)=bustt(choix1,choix2);
           taxi(index1,index2)=taxitt(choix1,choix2);
           train(index1,index2)=traintt(choix1,choix2);
       end
   end
   waitbar(index1/length(TZ2013));
end

%calcul de la distance ? vol d'oiseau d'une zone ? l'autre
distance_vol_oiseau=zeros(size(TZ2013));
for index=1:length(TZ2013),
    distance_vol_oiseau(:,index)=sqrt((X-X(index)).^2+(Y-Y(index)).^2)/1000;
end

%sauvegarde de tout ?a
save('transport_time','cars','bus','train','taxi', 'X','Y','distance_vol_oiseau','Area');

%% OPTIONNEL: calcul des densit;es d'emploi liss?es apr surface

emploi=zeros(size(TZ2013));
for index1=1:length(TZ2013),
    garde=(cars(index1,:)<=5);
    garde(cars(index1,:)==0)=0;
    garde(index1)=1;
    if sum(garde)>0,
       %pour avoir la densit? d'emploi
       emploi(index1)=sum(job_total(garde).*Area(garde))./sum(Area(garde));
       %pour avoir la somme des smplois ? moins d'un certain temps en
       %voiture
       %emploi(index1)=sum(job_total(garde));
    end
end

%% (pour tracer une carte)
trace_points_carte( cdata, X/1000,Y/1000,emploi);

