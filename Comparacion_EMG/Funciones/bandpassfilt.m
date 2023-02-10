% Función de filtro pasa bandas
% Esta función permite que se diseñe un filtro pasa bandas, 
% se requieren varios parametros, entre ellos:
% * fc1 = frecuencia de corte para el pasa altas;
% * fc2 frecuencia de corte para el pasa bajas;
% * fs = frecuencia de muestreo a la que se registraron los datos;
% * order = orden del filtro buscado;
% *data = datos que se busca filtrar

function data_filter = bandpassfilt(fc1, fc2, fs, order, data)
    [b,a] = butter(order,fc1/(fs/2),"high");
    data_1 = filter(b, a, data);

    [d,c] = butter(order,fc2/(fs/2),"low");
    data_2 = filter(d, c, data_1);

    data_filter = data_2;
end