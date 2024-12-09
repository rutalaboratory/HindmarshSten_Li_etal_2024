function relativeAngle = computeRelativeAngle(selfPos, otherPos, facingDir)

    % Compute the normalized position of other fly relative to self
    normSelfPos = otherPos-selfPos;
    
    % Compute the absolute angle of the other fly from self
    absAngle = atan2(normSelfPos(2, :), normSelfPos(1, :));
    
    % Compute the facing angle of other fly from self
    relativeAngle = circ_dist(absAngle, facingDir);
    
end