function closestValue = findClosestValue(v, a)

    [~, closestIndex] = min(abs(a - v));
    closestValue = v(closestIndex);

end