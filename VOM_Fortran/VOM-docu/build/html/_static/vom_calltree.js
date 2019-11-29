function vom(node, onCompleteCallback)
{
   var myobj = { label: "transpmodel", id: "transpmodel", href: "transpmodel.html", target:"basefrm" };
   var tmpNode2 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "sce", id: "sce", href: "sce.html", target:"basefrm" };
   var tmpNode32 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
