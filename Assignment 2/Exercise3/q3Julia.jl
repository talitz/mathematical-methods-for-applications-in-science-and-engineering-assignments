########################## Question 2 ############################

############################ PART A ##############################
function gram_schmidt(A)
    n = size(A,2); # assume A is m×n
    m = size(A,1);
    R = zeros(n,n);
    Q = zeros(m,n);
    R[1,1] = norm(A[:,1]);
    Q[:,1] = A[:,1]/R[1,1];

    for i=2:n
        z = A[:,i];
        for j=1:i-1
            column_q_j = Q[:,j]';
            row_a_i = A[:,i];
            R[j,i] = vecdot(column_q_j , row_a_i)
            z = z - R[j,i]*Q[:,j];
        end
        R[i,i] = norm(z);
        Q[:,i] = z/R[i,i];
    end

    return Q,R
end

function modified_gram_schmidt(A)
    n = size(A,2); # assume A is m×n
    m = size(A,1);
    R = zeros(n,n);
    Q = zeros(m,n);
    R[1,1] = norm(A[:,1]);
    Q[:,1] = A[:,1]/R[1,1];

    for i=2:n
        Q[:,i] = A[:,i];
        for j=1:i-1
            R[j,i] = vecdot(Q[:,j]',Q[:,i]);
            Q[:,i] = Q[:,i] - R[j,i]*Q[:,j];
        end
        R[i,i] = norm(Q[:,i]);
        Q[:,i] = Q[:,i]/R[i,i];
    end

    return Q,R
end

############################ PART B ##############################
A_1 = [1.0 1.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];
A_e = [1.0 1.0 1.0; 1e-10 0.0 0.0; 0.0 1e-10 0.0; 0.0 0.0 1e-10];

# QR factorization of the matrix A
Q1,R1 = gram_schmidt(A_1);
Q2,R2 = gram_schmidt(A_e);
Q3,R3 = modified_gram_schmidt(A_1);
Q4,R4 = modified_gram_schmidt(A_e);
println("              PART B              ")
println("QR factorization of the matrix A");
println("Q1 is: ", Q1)
println("R1 is: ", R1)
println("Q2 is: ", Q2)
println("R2 is: ", R2)
println("Q3 is: ", Q3)
println("R3 is: ", R3)
println("Q4 is: ", Q4)
println("R4 is: ", R4)
println("")
############################ PART C ##############################
get_id_matrix = P -> eye(size(P,2))
get_norm = A -> vecnorm(A'*A-get_id_matrix(A))

F1 = get_norm(Q1);
F2 = get_norm(Q2);
F3 = get_norm(Q3);
F4 = get_norm(Q4);
println("              PART C              ")
println("# Calculating Fᵢ = ||QᵀᵢQᵢ-I|| ")
println("Gram Schmidt")
println("F1 is: ", F1)
println("F2 is: ", F2)
println("Gram Schmidt Modified")
println("F3 is: ", F3)
println("F4 is: ", F4)
