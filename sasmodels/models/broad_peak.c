double Iq(double q, double porod_scale, double porod_exp, 
    double peak_scale, double correlation_length,
    double peak_pos, double width_exp, double shape_exp)
{
    double z;
    double inten;

    z = fabs(q - peak_pos) * correlation_length;
    if(q ==0 || z ==0){printf("division by zero");}
    else{
        inten = pow(porod_scale / pow(q, porod_exp)
                 + peak_scale / (1 + pow(z, width_exp)), shape_exp);
        return inten;
    }   
}