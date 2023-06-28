function [sig] = ref_constract(t, freq,phase, nb_Harmonic)
    for h_i = 1:nb_Harmonic
        sig(2*h_i-1,:) = sin(2*pi*h_i*freq*t+h_i*phase);
        sig(2*h_i,:) = cos(2*pi*h_i*freq*t+h_i*phase);
    end
end