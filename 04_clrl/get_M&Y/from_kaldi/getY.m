% labels = load('task1_label.txt');
% Y = onehot(labels);
% save Y Y

%% =================================================

function Y = getY(labelfilepath)
    labels = load(labelfilepath);
    Y = onehot(labels);
end

