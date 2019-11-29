function transpmodel_daily_step(node, onCompleteCallback)
{
   var myobj = { label: "vom_write_dayyear", id: "vom_write_dayyear", href: "vom_write_dayyear.html", target:"basefrm" };
   var tmpNode44 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_add_yearly", id: "vom_add_yearly", href: "vom_add_yearly.html", target:"basefrm" };
   var tmpNode46 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_save_dayyear", id: "vom_save_dayyear", href: "vom_save_dayyear.html", target:"basefrm" };
   var tmpNode47 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_adapt_foliage", id: "vom_adapt_foliage", href: "vom_adapt_foliage.html", target:"basefrm" };
   var tmpNode48 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_adapt_roots", id: "vom_adapt_roots", href: "vom_adapt_roots.html", target:"basefrm" };
   var tmpNode49 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
