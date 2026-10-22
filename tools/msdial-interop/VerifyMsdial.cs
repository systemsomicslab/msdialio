// Standalone .NET Framework probe using the user's existing MS-DIAL assemblies.
// Compile with PowerShell Add-Type -OutputType ConsoleApplication.
using System;
using System.IO;
using System.Reflection;
using System.Collections.Generic;
using System.Threading.Tasks;

class VerifyMsdial {
    static string assemblies;
    [STAThread]
    static int Main(string[] args) {
        if (args.Length < 3 || args.Length % 2 != 1) {
            Console.Error.WriteLine("Usage: VerifyMsdial ASSEMBLY_FOLDER FORMAT FILE [FORMAT FILE ...]");
            return 2;
        }
        assemblies = args[0];
        AppDomain.CurrentDomain.AssemblyResolve += delegate(object sender, ResolveEventArgs e) {
            string path = Path.Combine(assemblies, new AssemblyName(e.Name).Name + ".dll");
            return File.Exists(path) ? Assembly.LoadFrom(path) : null;
        };
        Assembly common = Assembly.LoadFrom(Path.Combine(assemblies, "Common.dll"));
        Assembly core = Assembly.LoadFrom(Path.Combine(assemblies, "MsdialCore.dll"));
        Type handler = common.GetType("CompMs.Common.MessagePack.MessagePackDefaultHandler", true);
        int failures = 0;
        for (int i = 1; i < args.Length; i += 2) {
            string kind = args[i], path = args[i + 1];
            try {
                if (kind == "dataset" || kind == "dataset-model") {
                    Assembly lcms = Assembly.LoadFrom(Path.Combine(assemblies, "MsdialLcMsApi.dll"));
                    Type storage = lcms.GetType("CompMs.MsdialLcMsApi.DataObj.MsdialLcmsDataStorage", true);
                    object serializer = storage.GetProperty("Serializer").GetValue(null, null);
                    string folder = Path.GetDirectoryName(Path.GetFullPath(path));
                    Type managerType = core.GetType("CompMs.MsdialCore.Parser.DirectoryTreeStreamManager", true);
                    using (IDisposable manager = (IDisposable)Activator.CreateInstance(managerType, new object[] {folder})) {
                        Task task = (Task)serializer.GetType().GetMethod("LoadAsync").Invoke(serializer, new object[] {manager, Path.GetFileName(path), folder, ""});
                        task.GetAwaiter().GetResult();
                        if (kind == "dataset-model") {
                            object data = task.GetType().GetProperty("Result").GetValue(task, null);
                            Assembly reactive = Assembly.LoadFrom(Path.Combine(assemblies, "ReactiveProperty.dll"));
                            Type brokerType = reactive.GetType("Reactive.Bindings.Notifiers.MessageBroker", true);
                            object broker = brokerType.GetField("Default").GetValue(null);
                            Assembly gui = Assembly.LoadFrom(Path.Combine(assemblies, "MSDIAL.exe"));
                            Type modelType = gui.GetType("CompMs.App.Msdial.Model.Core.DatasetModel", true);
                            Console.WriteLine("STAGE constructing DatasetModel");
                            object model = Activator.CreateInstance(modelType, new object[] {data, broker});
                            Console.WriteLine("PASS DatasetModel construction");
                        }
                    }
                    Console.WriteLine("PASS full dataset " + path);
                    continue;
                }
                Type type;
                string method = "LoadFromStream";
                if (kind == "mddata") {
                    Assembly lcms = Assembly.LoadFrom(Path.Combine(assemblies, "MsdialLcMsApi.dll"));
                    type = lcms.GetType("CompMs.MsdialLcMsApi.DataObj.MsdialLcmsDataStorage", true);
                } else if (kind == "msp2") {
                    type = common.GetType("CompMs.Common.Components.MoleculeMsReference", true);
                    method = "LoadLargerListFromStream";
                } else if (kind == "dbs-storage") {
                    type = core.GetType("CompMs.MsdialCore.DataObj.DataBaseStorage", true);
                } else if (kind == "arf2") {
                    type = core.GetType("CompMs.MsdialCore.DataObj.AlignmentResultContainer", true);
                } else if (kind == "pai2" || kind == "arf") {
                    string name = kind == "arf" ? "AlignmentChromPeakFeature" : "ChromatogramPeakFeature";
                    type = typeof(List<>).MakeGenericType(core.GetType("CompMs.MsdialCore.DataObj." + name, true));
                    if (kind == "arf") method = "LoadLargerListFromStream";
                } else {
                    throw new ArgumentException("Unknown format: " + kind);
                }
                using (Stream stream = File.OpenRead(path)) {
                    handler.GetMethod(method).MakeGenericMethod(type).Invoke(null, new object[] {stream});
                }
                Console.WriteLine("PASS " + path);
            } catch (Exception e) {
                while (e.InnerException != null) e = e.InnerException;
                Console.WriteLine("FAIL " + path + ": " + e);
                failures++;
            }
        }
        return failures == 0 ? 0 : 1;
    }
}
