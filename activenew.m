function [A,b,solution]=activSet(A,solution,b,diagM,ele_table)
    pen_par=100;
    lambda=A*solution-b;
    contact_force=lambda;
    contact_force=-contact_force.*diagM;
    active_set=[];
    dof_touched=zeros(total_dof,1);
    for e=1:No_Els
        map = ele_table(e,:);
        for dof_index=map
            if dof_touched(dof_index)==0
               dof_touched(dof_index)=1;
            else
                continue;
            end
            [x_g,y_g]=node2pos(dof_index);
            obstacle_value=gfunc(x_g,y_g);
            solution_value=solution(dof_index);
            active_cond=lambda(dof_index) + ...
            pen_par*diagM(dof_index) * ...
            (solution_value-obstacle_value)<0;
            if active_cond
                active_set=[active_set;dof_index];
                solution(dof_index) = obstacle_value;
                A(dof_index,:)=0;A(dof_index,dof_index)=1;
                b(dof_index)=obstacle_value;
                lambda(dof_index) = 0;
            end

        end
    end
end