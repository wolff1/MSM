function [ output_args ] = bounding_box_test(domain, q, d)
%bounding_box_test - determine size of largest bounding box
%                    such that no more than q:=sqrt(N) particles
%                    can be "contained" within the box.
    if domain == 0
        r = [0.0 0.0;
             0.0 1.0;
             1.0 0.0;
             1.0 1.0;
%              0.5 0.5
            ];
    else
        r = rand(q*q,d)*domain - domain/2.0;
    end
    N = length(r);
    q = floor(sqrt(N));
    d = size(r,2);

    LoCorner = min(r);
    HiCorner = max(r);
    amax = norm(HiCorner-LoCorner) + 1.0; % good idea / bad idea?

    % Calculate distance between each pair of particles
% O(N^2)
    D = zeros(N,N);
    for i = 1:N
        ri = r(i,1:d);
        for j = i+1:N
            rj = r(j,1:d);
            D(i,j) = norm(ri-rj);
            D(j,i) = D(i,j);
        end
    end

% ATTEMPT #3
    % Find q closest neighbors to particle i
% O(N^2 log(N))
    QNN = zeros(N,q);
    QNNI = zeros(N,q);
    for i = 1:N
        [v,idx] = sort(D(i,:));
        QNN(i,:) = v(2:q+1);
        QNNI(i,:) = idx(2:q+1);
    end

    % The distance from r(i) to x(q) defines the value of a for each r(i)
% O(N)
    a = amax;
    u = zeros(1,d);
    bin_origin = 0;
    for i = 1:N
        % Get positions of q nearest neighbors to particle i
        x = r(QNNI(i,:),1:d);

        v = x(q,1:d)-r(i,1:d);
        radius = sqrt(4.0/d)*norm(v);
        if radius < a
            a = radius;
            u = v;
            bin_origin = i;
        end
    end

    u = u / norm(u);                % Source Vector
    e = ones(1,d) / sqrt(d);        % Target Vector

    %   Build rotation matrix to rotate u to be aligned with e
    if d == 2
        theta = -acosd(u*e');
        R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    else
% FIXME - THIS NEEDS TO BE DOUBLE CHECKED!!!
%         Rx = [1.0 0.0 0.0; 0.0 cosd(theta) -sind(theta); 0.0 sind(theta) cosd(theta)];

        e(3) = 0.0;                 % [sqrt(d) sqrt(d) 0.0] -> xy-plane
        theta = -acosd(u*e');
        R = [cosd(theta) -sind(theta) 0.0; sind(theta) cosd(theta) 0.0; 0.0 0.0 1.0];

        e(3) = e(1);
        e(2) = 0.0;                 % [sqrt(d) 0.0 sqrt(d)] -> xz-plane
        theta = -acosd(u*e');
%O(d^3)
        R = R*[cosd(theta) 0.0 sind(theta); 0.0 1.0 0.0; -sind(theta) 0.0 cosd(theta)];
    end

    % Rotate particles about the axes so that bin-ing can proceed
% O(d^2N)
    for i = 1:N
        r(i,1:d) = R*(r(i,1:d)-r(bin_origin,1:d))';
    end

    % With a and new basis, set up bins
    LoCorner = min(r);
    HiCorner = max(r);
    LoIdx = -floor(LoCorner/a);
    HiIdx = ceil(HiCorner/a);
    bins = LoIdx + HiIdx + ones(1,d);
    bin = zeros(bins);

    % Place particles into bins
% O(N)
    if d == 2
        for i = 1:N
            bidx = floor(r(i,1:d)/a) + LoIdx + ones(1,d);
            bin(bidx(1),bidx(2)) = bin(bidx(1),bidx(2)) + 1;
        end
        max_in_bins = norm(reshape(bin,1,bins(1)*bins(2)),inf);
        sum_in_bins = sum(reshape(bin,1,bins(1)*bins(2)));
    else
        for i = 1:N
            bidx = floor(r(i,1:d)/a) + LoIdx + ones(1,d);
            bin(bidx(1),bidx(2),bidx(3)) = bin(bidx(1),bidx(2),bidx(3)) + 1;
        end
        max_in_bins = norm(reshape(bin,1,bins(1)*bins(2)*bins(3)),inf);
        sum_in_bins = sum(reshape(bin,1,bins(1)*bins(2)*bins(3)));
    end
    fprintf('amax = %f, a=%f, N=%d, q=%d, max-in-bin=%d, N''=%d, a pct = %f\n', amax-1.0, a, N, q, max_in_bins, sum_in_bins, a*100.0/(amax-1.0));

%-----------------------------------------------------------------------
% % ATTEMPT #1
%     % Find q closest neighbors to particle i, O(N^(3/2) log(N^(1/2)))
%     QNN = zeros(N,q+1);
%     QNNI = zeros(N,q+1);
%     for i = 1:N
%         [v,idx] = sort(D(i,:));
%         QNN(i,:) = v(1:q+1);
%         QNNI(i,:) = idx(1:q+1);
%     end
% %    QNN
% %    QNNI
% 
%     % Minimum side length must be at least distance to nearest neighbor
%     % In 2D, "in box" is solution of 4 linear inequalities (4 lines)
%     %   "between" 4 points (corners of the box)
%     % In 3D, "in box" is solution of 6 linear inequalities (6 planes)
%     %   "between" 8 points (corners of the cube)
%     % Need to determine the corners of box (2D) / cube (3D) ?
%     %   Probably in order to form inequalities
%     % Particles lying on "lower edges" are in (?)
%     % Particles lying on "upper edges" are out (?)
% 
% %   X         X
% %        X     
% %   X         X
% 
% %   ia <= x < (i+1)a -> in box i
% 
%     % MAX(QNN(:,q))  -> distance b/t particles s.t. box contains q
%     % MIN(QNN(:,q+1))-> distance b/t particles s.t. box contains q
% 
%     a = amax;
%     D2 = zeros(q+1,1);
%     for i = 1:N
%         % Get positions of q nearest neighbors to particle i
%         x = r(QNNI(i,:),1:d);
%         for j = 1:q+1
%             D2(j) = norm(x(j,:)-LoCorner);
%         end
%         [~,idx2] = sort(D2);
% %        x(idx2,:) % particle positions "in order" from bottom corner to furtherst
% 
%         %min_coord = min(min(x(idx2,:)));
%         %max_coord = max(max(x(idx2,:)));
%         %radius = max_coord-min_coord;
%         radius = norm(x(idx2(q+1),:) - x(idx2(1),:));
%         a = min(a, radius);
%     end
%     % For each particle, find *largest* box that contains it and q-1 nearest
%     % neighbors
%     % Take the smallest such box.
%     a = max(1.0,a); % 1.0 <= a <= amax
%     fprintf('amax = %f, a=%f, N=%d, q=%d\n', amax, a, N, q);

% % ATTEMPT #2
%     % Find q-1 closest neighbors to particle i, O(N^(3/2) log(N^(1/2)))
%     QNN = zeros(N,q);
%     QNNI = zeros(N,q);
%     for i = 1:N
%         [v,idx] = sort(D(i,:));
%         QNN(i,:) = v(1:q);
%         QNNI(i,:) = idx(1:q);
%     end
% %    QNN
% %    QNNI
% 
%     % find center of particle i and q-1 nearest neighbors
%     % find distance from center to nearest particle outside of q cloud
%     a = amax;
%     for i = 1:N
%         % Get positions of q-1 nearest neighbors to particle i
%         x = r(QNNI(i,:),1:d);
% 
%         % Find center of q particles which includes particle i
%         c = sum(x)/q;
% %         c
% 
% %         % Find radius from center to furtherst particle in x
% %         radius = 0.0;
% %         for j = 1:q
% %             radius = max(radius, norm(c-x(j,:)));
% %         end
% %         radius
% 
%         % Find distance from center to q+1st nearest particle
%         D3 = zeros(N,1);
%         for j = 1:N
%             D3(j) = norm(c-r(j,:));
%         end
% %         D3
%         [~,idx3] = sort(D3);
%         qp1 = r(idx3(q+1),:);
% 
%         % distance from c to qp1 is half-diagonal of bounding box
%         radius = norm(qp1-c);
% 
%         sqrt(2.0)*radius
%         a = min(a, sqrt(2.0)*radius);
%     end
%     a = max(1.0, a);
%     a
end
