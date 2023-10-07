function Y = getY(labelfilepath)
    labels = load(labelfilepath);
    Y = onehot(labels);
end

