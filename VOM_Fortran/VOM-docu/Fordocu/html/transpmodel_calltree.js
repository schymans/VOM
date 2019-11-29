function transpmodel(node, onCompleteCallback)
{
   var myobj = { label: "transpmodel_init", id: "transpmodel_init", href: "transpmodel_init.html", target:"basefrm" };
   var tmpNode22 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_daily_init", id: "vom_daily_init", href: "vom_daily_init.html", target:"basefrm" };
   var tmpNode25 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_hourly_init", id: "vom_hourly_init", href: "vom_hourly_init.html", target:"basefrm" };
   var tmpNode26 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_gstom", id: "vom_gstom", href: "vom_gstom.html", target:"basefrm" };
   var tmpNode27 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_subhourly_init", id: "vom_subhourly_init", href: "vom_subhourly_init.html", target:"basefrm" };
   var tmpNode28 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode28.isLeaf = true; 
   var myobj = { label: "vom_rootuptake", id: "vom_rootuptake", href: "vom_rootuptake.html", target:"basefrm" };
   var tmpNode29 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_mqss", id: "vom_mqss", href: "vom_mqss.html", target:"basefrm" };
   var tmpNode30 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode30.isLeaf = true; 
   var myobj = { label: "vom_tissue_water_et", id: "vom_tissue_water_et", href: "vom_tissue_water_et.html", target:"basefrm" };
   var tmpNode31 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_subhourly", id: "vom_subhourly", href: "vom_subhourly.html", target:"basefrm" };
   var tmpNode33 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_add_hourly", id: "vom_add_hourly", href: "vom_add_hourly.html", target:"basefrm" };
   var tmpNode39 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_add_daily", id: "vom_add_daily", href: "vom_add_daily.html", target:"basefrm" };
   var tmpNode40 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_write_hourly", id: "vom_write_hourly", href: "vom_write_hourly.html", target:"basefrm" };
   var tmpNode41 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_check_water", id: "vom_check_water", href: "vom_check_water.html", target:"basefrm" };
   var tmpNode42 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode42.isLeaf = true; 
   var myobj = { label: "transpmodel_daily_step", id: "transpmodel_daily_step", href: "transpmodel_daily_step.html", target:"basefrm" };
   var tmpNode43 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "transpmodel_last_step", id: "transpmodel_last_step", href: "transpmodel_last_step.html", target:"basefrm" };
   var tmpNode50 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode50.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
