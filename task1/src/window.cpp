#include <window.hpp>
#include <solution.hpp>
#include <qlayout.h>
#include <qlabel.h>
#include <qgroupbox.h>
#include <qspinbox.h>
#include <qpushbutton.h>
#include <qchart.h>
#include <qchartview.h>
#include <qvalueaxis.h>
#include <qlineseries.h>
#include <qlegendmarker.h>

constexpr double x0def{4};
constexpr double t0def{0};
constexpr double tkdef{10};
constexpr double dtdef{0.05};
constexpr std::size_t maxiterdef{100};

using namespace QtCharts;

auto make_line_series(QString const &name, QColor const &color)
{
    auto series = new QLineSeries;
    series->setName(name);
    series->setPen(QPen(color, 3));
    return series;
}

auto make_axis(QString const &name)
{
    auto axis = new QValueAxis;
    axis->setGridLineVisible(true);
    axis->setTitleText(name);
    axis->setLabelFormat("%.1f");
    axis->setMinorTickCount(0.5);
    return axis;
}

auto make_chartview(QAbstractAxis *xaxis, QAbstractAxis *yaxis, QAbstractSeries *series)
{
    auto chart = new QChart;
    chart->addAxis(xaxis, Qt::AlignmentFlag::AlignBottom);
    chart->addAxis(yaxis, Qt::AlignmentFlag::AlignLeft);
    chart->addSeries(series);
    chart->legend()->markers(series).front()->setVisible(false);
    series->attachAxis(xaxis);
    series->attachAxis(yaxis);
    auto view = new QChartView(chart);
    view->setRenderHint(QPainter::RenderHint::Antialiasing);
    return view;
}

template <double node::*field>
class node_iterator
{
    node *_ptr;

public:
    node_iterator(node *ptr) : _ptr{ptr} {}
    node_iterator(node_iterator const &) = default;
    double operator*() const { return _ptr->*field; }
    bool operator!=(node_iterator const &other) const { return _ptr != other._ptr; }
    bool operator==(node_iterator const &other) const { return _ptr == other._ptr; }
    node_iterator &operator++()
    {
        ++_ptr;
        return *this;
    }
    node_iterator operator++(int)
    {
        auto iter = *this;
        ++_ptr;
        return iter;
    }

public:
    using value_type = double;
};

class data_drawer : public QObject
{
    QLineSeries *_series;
    QValueAxis *_xaxis, *_yaxis;
    QChartView *_view;

public:
    data_drawer(QString const &xname, QString const &yname, QColor color, QObject *parent) : QObject(parent)
    {
        _series = make_line_series(yname, color);
        _xaxis = make_axis(xname);
        _yaxis = make_axis(yname);
        _view = make_chartview(_xaxis, _yaxis, _series);
        _view->setMinimumHeight(200);
    }
    QChartView *get_view() const { return _view; }
    template <double node::*f>
    void draw(node_iterator<&node::t> xiter, node_iterator<f> yiter, std::size_t count, double xmin, double xmax, double ymin, double ymax)
    {
        _series->clear();
        for (; count > 0; --count)
        {
            _series->append(*xiter++, *yiter++);
        }
        double dy = (ymax - ymin) * 1e-2;
        _xaxis->setRange(xmin, xmax);
        _yaxis->setRange(ymin - dy, ymax + dy);
        _view->update();
    }
};

class appmodel : public QObject
{
    solution _sol;
    data_drawer *_xdrawer, *_pdrawer, *_udrawer;
    QLabel *_lbl;
    double _x0{x0def};
    double _t0{t0def};
    double _tk{tkdef};
    double _dt{dtdef};
    double _eps{1e-3};
    std::size_t _maxiter{maxiterdef};
    char const *_filepath{"log.txt"};

public:
    appmodel(data_drawer *xdrawer, data_drawer *pdrawer, data_drawer *udrawer, QLabel *lbl, QObject *parent)
        : QObject(parent), _xdrawer{xdrawer}, _pdrawer{pdrawer}, _udrawer{udrawer}, _lbl{lbl}
    {
    }
    void set_start_coordinate(double x0)
    {
        _x0 = x0;
    }
    void set_start_time(double t0)
    {
        _t0 = t0;
    }
    void set_time_interval(double interval)
    {
        _tk = _t0 + interval;
    }
    void set_delta_time(double dt)
    {
        _dt = dt;
    }
    void set_iterations_count(int count)
    {
        _maxiter = static_cast<std::size_t>(count);
    }
    void set_precision(double eps)
    {
        _eps = eps;
    }
    void compute(bool)
    {
        _sol = solve(_x0, _t0, _tk, _dt, _eps, _maxiter, _filepath);
        _lbl->setText("Значение ф-ла " + QString::number(_sol.funcval, 'f', 3));
        double tmin = _sol.nodes.front().t, tmax = _sol.nodes.back().t;
        std::size_t size = _sol.nodes.size();
        auto begin_ptr = _sol.nodes.data(), end_ptr = begin_ptr + size;
        auto xiter_pair = std::minmax_element(node_iterator<&node::x>{begin_ptr}, node_iterator<&node::x>{end_ptr});
        auto piter_pair = std::minmax_element(node_iterator<&node::p>{begin_ptr}, node_iterator<&node::p>{end_ptr});
        _xdrawer->draw(node_iterator<&node::t>{begin_ptr}, node_iterator<&node::x>{begin_ptr}, size, tmin, tmax, *xiter_pair.first, *xiter_pair.second);
        _pdrawer->draw(node_iterator<&node::t>{begin_ptr}, node_iterator<&node::p>{begin_ptr}, size, tmin, tmax, *piter_pair.first, *piter_pair.second);
        _udrawer->draw(node_iterator<&node::t>{begin_ptr}, node_iterator<&node::u>{begin_ptr}, size, tmin, tmax, -1, 1);
    }
};

constexpr int btnwidth{100};
constexpr int btnheight{30};
constexpr Qt::Alignment alignment{Qt::AlignmentFlag::AlignRight + Qt::AlignmentFlag::AlignVCenter};

auto make_label(QString const &text)
{
    auto lbl = new QLabel(text);
    lbl->setFixedHeight(btnheight);
    return lbl;
}

appwindow::appwindow()
{
    auto xdrawer = new data_drawer("t", "x", QColor(Qt::GlobalColor::blue), this);
    auto pdrawer = new data_drawer("t", "p", QColor(Qt::GlobalColor::green), this);
    auto udrawer = new data_drawer("t", "u", QColor(Qt::GlobalColor::red), this);
    auto lbl = make_label("Значение ф-ла ---");
    auto model = new appmodel(xdrawer, pdrawer, udrawer, lbl, this);
    setWindowTitle("Окно приложения");
    setCentralWidget(new QWidget);
    auto grid = new QGridLayout;
    grid->setRowStretch(0, 0);
    grid->setRowStretch(1, 0);
    grid->setRowStretch(2, 1);
    grid->setRowStretch(3, 1);
    grid->setRowStretch(4, 1);
    centralWidget()->setLayout(grid);

    auto spb_t0 = new QDoubleSpinBox;
    spb_t0->setValue(t0def);
    spb_t0->setSingleStep(1);
    spb_t0->setMinimum(0);
    spb_t0->setFixedSize(btnwidth, btnheight);
    auto spb_interval = new QDoubleSpinBox;
    spb_interval->setValue(tkdef - t0def);
    spb_interval->setSingleStep(1);
    spb_interval->setMinimum(1);
    spb_interval->setFixedSize(btnwidth, btnheight);
    auto spb_dt = new QDoubleSpinBox;
    spb_dt->setValue(dtdef);
    spb_dt->setSingleStep(1e-1);
    spb_dt->setMinimum(1e-2);
    spb_dt->setMaximum(spb_interval->value());
    spb_dt->setFixedSize(btnwidth, btnheight);
    auto spb_x0 = new QDoubleSpinBox;
    spb_x0->setMinimum(-1000);
    spb_x0->setMaximum(1000);
    spb_x0->setValue(4);
    spb_x0->setSingleStep(1);
    spb_x0->setFixedSize(btnwidth, btnheight);
    auto spb_maxiter = new QSpinBox;
    spb_maxiter->setMinimum(1);
    spb_maxiter->setMaximum(1000);
    spb_maxiter->setValue(maxiterdef);
    spb_maxiter->setFixedSize(btnwidth, btnheight);
    auto btn = new QPushButton("Рассчитать");
    btn->setFixedSize(btnwidth, btnheight);

    connect(btn, &QPushButton::clicked, model, &appmodel::compute);
    connect(spb_t0, QOverload<double>::of(&QDoubleSpinBox::valueChanged), model, &appmodel::set_start_time);
    connect(spb_interval, QOverload<double>::of(&QDoubleSpinBox::valueChanged), model, &appmodel::set_time_interval);
    connect(spb_interval, QOverload<double>::of(&QDoubleSpinBox::valueChanged), spb_dt, &QDoubleSpinBox::setMaximum);
    connect(spb_dt, QOverload<double>::of(&QDoubleSpinBox::valueChanged), model, &appmodel::set_delta_time);
    connect(spb_x0, QOverload<double>::of(&QDoubleSpinBox::valueChanged), model, &appmodel::set_start_coordinate);
    connect(spb_maxiter, QOverload<int>::of(&QSpinBox::valueChanged), model, &appmodel::set_iterations_count);

    grid->addWidget(make_label("Нач. время"), 0, 0, alignment);
    grid->addWidget(spb_t0, 0, 1);
    grid->addWidget(make_label("Интервал вр."), 0, 2, alignment);
    grid->addWidget(spb_interval, 0, 3);
    grid->addWidget(make_label("Шаг интегр."), 0, 4, alignment);
    grid->addWidget(spb_dt, 0, 5);
    grid->addWidget(make_label("Нач. коорд."), 0, 6, alignment);
    grid->addWidget(spb_x0, 0, 7);
    grid->addWidget(make_label("Кол-во итер."), 0, 8, alignment);
    grid->addWidget(spb_maxiter, 0, 9);
    grid->addWidget(lbl, 1, 0, 1, 2);
    grid->addWidget(btn, 1, 2, 1, 8, Qt::AlignmentFlag::AlignCenter);
    grid->addWidget(xdrawer->get_view(), 2, 0, 1, 10);
    grid->addWidget(pdrawer->get_view(), 3, 0, 1, 10);
    grid->addWidget(udrawer->get_view(), 4, 0, 1, 10);
}