
void generate_gaussian()
{
    RooWorkspace ws("ws");
    ws.factory("Gaussian::gaus(x[-10,10], mean[0], sigma[1])");
    ws.factory("PROD::model(gaus)");

    RooDataSet* data = ws.pdf("model")->generate(*ws.var("x"), 1000);
    RooPlot *frame = ws.var("x")->frame();
    data           ->plotOn(frame);
    ws.pdf("model")->plotOn(frame);
    frame->Draw();
}
