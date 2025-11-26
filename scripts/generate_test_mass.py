import ROOT as rt

input_file = rt.TFile("../input/lhc16q_final.root")
input_list = input_file.Get("h-k0")

input_dist = input_list.FindObject("fK0Dist")

input_dist.GetAxis(0).SetRangeUser(1.0, 2.999)

mass_dist = input_dist.Projection(3).Clone("test_mass_dist")

output_file = rt.TFile("test_mass.root", "RECREATE")
output_file.cd()
mass_dist.Write()


