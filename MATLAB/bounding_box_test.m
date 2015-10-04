function [ output_args ] = bounding_box_test(domain, q, d)
%bounding_box_test - determine size of largest bounding box
%                    such that no more than q:=sqrt(N) particles
%                    can be "contained" within the box.
%     r = [0.0 0.0;
%          0.0 1.0;
%          1.0 0.0;
%          1.0 1.0;
% %         0.5 0.5
%         ];
    r = rand(q*q,d)*domain - domain/2.0;
    r
    N = length(r);
    q = floor(sqrt(N));
    d = size(r,2);

    LoCorner = min(r);
    HiCorner = max(r);
    amax = norm(HiCorner-LoCorner);

    % Calculate distance between each pair of particles, O(N^2)
    D = zeros(N,N);
    for i = 1:N
        ri = r(i,1:d);
        for j = i+1:N
            rj = r(j,1:d);
            D(i,j) = norm(ri-rj);
            D(j,i) = norm(ri-rj);
        end
    end
%    D

    % Find q closest neighbors to particle i, O(N^(3/2) log(N^(1/2)))
    QNN = zeros(N,q+1);
    QNNI = zeros(N,q+1);
    for i = 1:N
        [v,idx] = sort(D(i,:));
        QNN(i,:) = v(1:q+1);
        QNNI(i,:) = idx(1:q+1);
    end
%    QNN
%    QNNI

    % Minimum side length must be at least distance to nearest neighbor
    % In 2D, "in box" is solution of 4 linear inequalities (4 lines)
    %   "between" 4 points (corners of the box)
    % In 3D, "in box" is solution of 6 linear inequalities (6 planes)
    %   "between" 8 points (corners of the cube)
    % Need to determine the corners of box (2D) / cube (3D) ?
    %   Probably in order to form inequalities
    % Particles lying on "lower edges" are in (?)
    % Particles lying on "upper edges" are out (?)

%   X         X
%        X     
%   X         X

%   ia <= x < (i+1)a -> in box i

    % MAX(QNN(:,q))  -> distance b/t particles s.t. box contains q
    % MIN(QNN(:,q+1))-> distance b/t particles s.t. box contains q

    a = amax;
    D2 = zeros(q+1,1);
    for i = 1:N
        % Get positions of q nearest neighbors to particle i
        x = r(QNNI(i,:),1:d);
        for j = 1:q+1
            D2(j) = norm(x(j,:)-LoCorner);
        end
        [~,idx2] = sort(D2);
%        x(idx2,:) % particle positions "in order" from bottom corner to furtherst

        %min_coord = min(min(x(idx2,:)));
        %max_coord = max(max(x(idx2,:)));
        %radius = max_coord-min_coord;
        radius = norm(x(idx2(q+1),:) - x(idx2(1),:));
        a = min(a, radius);
    end
    % For each particle, find *largest* box that contains it and q-1 nearest
    % neighbors
    % Take the smallest such box.
    a = max(1.0,a); % 1.0 <= a <= amax
    fprintf('amax = %f, a=%f, N=%d, q=%d\n', amax, a, N, q);

    % find center of particle i and q-1 nearest neighbors
    % find distance from center to nearest particle outside of q cloud
end
