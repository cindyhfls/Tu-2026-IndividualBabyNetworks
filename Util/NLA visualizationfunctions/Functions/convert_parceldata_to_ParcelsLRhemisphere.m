function Parcels = convert_parceldata_to_ParcelsLRhemisphere(parceldata)
% parceldata is 59412 x 1 vector for cortical parcel assignment in fsLR-32k
% space
% Parcels is a struture with CtxL and CtxR and medial wall assigned 0.

Linds=with_without_mw_conversion('Lindfull');
Rinds=with_without_mw_conversion('Rindfull');
   
n_verts_per_hem = 32492;
Parcels.CtxL = zeros(n_verts_per_hem,1);% total vertices
Parcels.CtxR = zeros(n_verts_per_hem,1);
Parcels.CtxL(Linds) = parceldata(1:length(Linds));
Parcels.CtxR(Rinds-n_verts_per_hem) = parceldata(length(Linds)+1:end);

end